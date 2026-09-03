use anyhow::{bail, Context, Result};
use bio::io::fastq::{Reader as FastqReader, Record, Writer as FastqWriter};
use clap::{ArgAction, Parser};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use regex::Regex;
use std::collections::{BTreeMap, HashMap};
use std::fs::{create_dir_all, File};
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;
use bio::io::fastq::FastqRead;

/// CLI definition
#[derive(Parser, Debug)]
#[command(name = "fastq_barcode_splitter", author, version, about = "Split/combine FASTQ reads by barcodes per pool")]
struct Cli {
    /// Input directory containing FASTQ.gz files
    #[arg(short = 'i', long = "input")]
    input: PathBuf,

    /// Output directory (created if missing)
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Max number of pools processed in parallel (up to 4)
    #[arg(short = 't', long = "threads", default_value_t = 4)]
    threads: usize,

    /// Verbose warnings (prints every ID mismatch)
    #[arg(long = "verbose", action = ArgAction::SetTrue)]
    verbose: bool,

    /// Diagnose barcode frequencies (prints top 20 four-mers at [2..6] for R1 and R2 per pool)
    #[arg(long = "diagnose", action = ArgAction::SetTrue)]
    diagnose: bool,

    /// Reverse-complement the 4-mer slice from R2 (P5) before matching
    #[arg(long = "rc-p5", action = ArgAction::SetTrue)]
    rc_p5: bool,

    /// Reverse-complement the 4-mer slice from R1 (P7) before matching
    #[arg(long = "rc-p7", action = ArgAction::SetTrue)]
    rc_p7: bool,
}

/// Sample metadata: names and exact barcode pairs
#[derive(Clone)]
struct SampleDef {
    name: &'static str,
    p5: &'static [u8], // reverse read barcode (R2)
    p7: &'static [u8], // forward read barcode (R1)
}

/// Output stats per pool
#[derive(Default)]
struct PoolStats {
    pool: String,
    total_pairs: u64,
    unassigned: u64,
    per_sample_assigned: BTreeMap<String, u64>,
}

/// A pair of lane files (R1,R2)
#[derive(Clone)]
struct LanePair {
    r1: PathBuf,
    r2: PathBuf,
    lane: String,
}

/// Input set grouped by pool
#[derive(Default)]
struct PoolInput {
    pool_id: String,     // "03", "10", "11", "12"
    lanes: Vec<LanePair>,// (R1,R2) pairs per lane
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Sanity checks
    if !cli.input.is_dir() {
        bail!("--input must be a directory");
    }
    create_dir_all(&cli.output).context("creating output directory")?;

    // Discover input files and group per pool
    let pools = discover_inputs(&cli.input)?;
    if pools.is_empty() {
        bail!("No matching FASTQ.gz files found in {:?}", cli.input);
    }

    // Cap threads by number of pools and 4 as requested
    let num_jobs = pools.len().min(cli.threads).min(4);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_jobs)
        .build()
        .context("building rayon thread pool")?;

    // Prepare barcode dictionaries
    let control = SampleDef {
        name: "Control",
        p5: b"ATGC",
        p7: b"TCGC",
    };
    let experimental = vec![
        SampleDef { name: "WT_SpCas9", p5: b"ATGC", p7: b"AGTA" },
        SampleDef { name: "SpG",       p5: b"ATGC", p7: b"ACTT" },
        SampleDef { name: "SpRY",      p5: b"GCAA", p7: b"TCCA" },
    ];

    let exp_lookup: HashMap<(Vec<u8>, Vec<u8>), &'static str> = experimental
        .iter()
        .map(|s| ((s.p5.to_vec(), s.p7.to_vec()), s.name))
        .collect();

    // Process pools in parallel
    let results: Vec<PoolStats> = pool.install(|| {
        pools
            .par_iter()
            .map(|p| {
                process_pool(
                    p,
                    &cli.output,
                    &control,
                    &experimental,
                    &exp_lookup,
                    cli.verbose,
                    cli.diagnose,
                    cli.rc_p5,
                    cli.rc_p7,
                )
            })
            .collect()
    });

    // Print summary
    println!("\n===== Summary =====");
    for stats in results {
        println!("Pool {}:", stats.pool);
        println!("  Total read pairs processed: {}", stats.total_pairs);
        for (sample, count) in stats.per_sample_assigned.iter() {
            println!("    {}: {}/{}", sample, count, stats.total_pairs);
        }
        println!("  Unassigned: {}", stats.unassigned);
        println!();
    }

    Ok(())
}

/// Discover FASTQ.gz files and group by pool, pairing R1/R2 per lane
fn discover_inputs(input_dir: &Path) -> Result<Vec<PoolInput>> {
    // Example filenames:
    // expRW086_pool_03_S3_L001_1.fastq.gz
    // expRW086_pool_10_S10_L004_2.fastq.gz
    let re = Regex::new(r"^expRW086_pool_(\d{2})_S\d+_L00([1-4])_(1|2)\.fastq\.gz$")?;

    // pool_id -> lane_idx -> {r1,r2,lane_label}
    let mut map: HashMap<String, HashMap<String, (Option<PathBuf>, Option<PathBuf>, String)>> = HashMap::new();

    for entry in WalkDir::new(input_dir).into_iter().filter_map(|e| e.ok()) {
        if !entry.file_type().is_file() {
            continue;
        }

        let fname = entry.file_name().to_string_lossy().to_string();
        if let Some(caps) = re.captures(&fname) {
            let pool_id = caps.get(1).unwrap().as_str().to_string();
            let lane_idx = caps.get(2).unwrap().as_str().to_string();
            let read = caps.get(3).unwrap().as_str(); // "1" or "2"
            let lane_label = format!("L00{}", lane_idx);

            let pool_entry = map.entry(pool_id.clone()).or_default();
            let lane_entry = pool_entry
                .entry(lane_idx.clone())
                .or_insert((None, None, lane_label.clone()));

            let path = entry.path().to_path_buf();
            match read {
                "1" => lane_entry.0 = Some(path),
                "2" => lane_entry.1 = Some(path),
                _ => {}
            }
        }
    }

    // Build PoolInput vec, skipping lanes that don't have both R1 and R2
    let mut pools: Vec<PoolInput> = Vec::new();
    for (pool_id, lanes_map) in map.into_iter() {
        let mut pi = PoolInput { pool_id: pool_id.clone(), lanes: Vec::new() };
        for (lane_idx, (r1_opt, r2_opt, lane_label)) in lanes_map.into_iter() {
            match (r1_opt, r2_opt) {
                (Some(r1), Some(r2)) => {
                    pi.lanes.push(LanePair { r1, r2, lane: lane_label });
                }
                _ => {
                    eprintln!(
                        "[WARN] Skipping pool {} lane {}: missing R1 or R2 file",
                        pool_id, lane_idx
                    );
                }
            }
        }
        // Sort lanes by lane label (L001..L004) for deterministic processing
        pi.lanes.sort_by(|a, b| a.lane.cmp(&b.lane));
        if !pi.lanes.is_empty() {
            pools.push(pi);
        }
    }

    // Sort pools by numeric id
    pools.sort_by(|a, b| a.pool_id.cmp(&b.pool_id));
    Ok(pools)
}

/// Process a single pool: combine lanes, split by barcodes, write outputs
fn process_pool(
    pool: &PoolInput,
    out_dir: &Path,
    control: &SampleDef,
    experimental: &[SampleDef],
    exp_lookup: &HashMap<(Vec<u8>, Vec<u8>), &'static str>,
    verbose: bool,
    do_diag: bool,
    rc_p5: bool,
    rc_p7: bool,
) -> PoolStats {
    let mut stats = PoolStats {
        pool: pool.pool_id.clone(),
        ..Default::default()
    };

    // For diagnostics: frequency maps of observed 4-mers at [2..6]
    let mut freq_r1: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut freq_r2: HashMap<Vec<u8>, u64> = HashMap::new();

    // Prepare output writers
    let mut writers: HashMap<String, FastqWriter<GzEncoder<BufWriter<File>>>> = HashMap::new();
    let is_control_pool = pool.pool_id == "03";

    if is_control_pool {
        let out_path = out_dir.join("expRW086_control_1.fastq.gz");
        let w = make_writer(&out_path);
        writers.insert(control.name.to_string(), w);
    } else {
        for s in experimental {
            let fname = format!("expRW086_{}_pool_{}_1.fastq.gz", s.name, pool.pool_id);
            let out_path = out_dir.join(fname);
            let w = make_writer(&out_path);
            writers.insert(s.name.to_string(), w);
        }
    }

    // Process each lane pair
    for lane in &pool.lanes {
        let mut r1 = make_reader(&lane.r1);
        let mut r2 = make_reader(&lane.r2);

        let mut rec1 = Record::new();
        let mut rec2 = Record::new();

        // NEW READ LOOP USING rust-bio 4.x EOF PATTERN
        loop {
            // Read next records; rust-bio returns Ok(()) and we check is_empty() for EOF
            if let Err(e) = r1.read(&mut rec1) {
                eprintln!("[ERROR] reading R1 {}: {}", lane.r1.display(), e);
                break;
            }
            if let Err(e) = r2.read(&mut rec2) {
                eprintln!("[ERROR] reading R2 {}: {}", lane.r2.display(), e);
                break;
            }

            // EOF handling per rust-bio idiom
            let eof1 = rec1.is_empty();
            let eof2 = rec2.is_empty();
            if eof1 && eof2 {
                break; // clean EOF on both
            }
            if eof1 ^ eof2 {
                eprintln!(
                    "[WARN] Pool {} lane {}: R1/R2 length mismatch (EOF on one file). Stopping lane.",
                    pool.pool_id, lane.lane
                );
                break;
            }

            // From here, both rec1 and rec2 contain a pair
            stats.total_pairs += 1;

            // ID check (remove /1 or /2; use first token)
            let id1 = canonical_id(rec1.id());
            let id2 = canonical_id(rec2.id());
            if id1 != id2 {
                if verbose {
                    eprintln!(
                        "[WARN] Pool {} lane {}: R1/R2 ID mismatch -> R1:{} R2:{} (files: {}, {})",
                        pool.pool_id, lane.lane,
                        rec1.id(), rec2.id(),
                        lane.r1.display(), lane.r2.display()
                    );
                }
                stats.unassigned += 1;
                continue;
            }

            // Ensure we can slice [2..6] in both reads
            let seq1 = rec1.seq();
            let seq2 = rec2.seq();
            if seq1.len() < 6 || seq2.len() < 6 {
                stats.unassigned += 1;
                continue;
            }

            // Extract 4-mers with corrected mapping: R1 => P5, R2 => P7
            let mut p5 = upper4(&seq1[2..6]); // R1 carries P5
            let mut p7 = upper4(&seq2[2..6]); // R2 carries P7
            if rc_p5 {
                p5 = revcomp4(&p5);
            }
            if rc_p7 {
                p7 = revcomp4(&p7);
            }

            // Collect diagnostic frequencies if requested
            if do_diag {
                bump(&mut freq_r1, p7.clone());
                bump(&mut freq_r2, p5.clone());
            }

            if is_control_pool {
                // Control requires exact match of both barcodes
                if p7 == control.p7 && p5 == control.p5 {
                    if let Some(w) = writers.get_mut(control.name) {
                        if let Err(e) = w.write_record(&rec1) {
                            eprintln!("[ERROR] write failed: {}", e);
                        } else {
                            *stats.per_sample_assigned.entry(control.name.to_string()).or_insert(0) += 1;
                        }
                    }
                } else {
                    stats.unassigned += 1;
                }
            } else {
                // Experimental: exact match against any sample pair
                if let Some(sample_name) = exp_lookup.get(&(p5.to_vec(), p7.to_vec())) {
                    if let Some(w) = writers.get_mut(&sample_name.to_string()) {
                        if let Err(e) = w.write_record(&rec1) {
                            eprintln!("[ERROR] write failed: {}", e);
                        } else {
                            *stats.per_sample_assigned.entry(sample_name.to_string()).or_insert(0) += 1;
                        }
                    }
                } else {
                    stats.unassigned += 1;
                }
            }
        }
    }

    // Print diagnostic frequencies (Top-20) if requested
    if do_diag {
        println!("\n[Diag] Pool {} barcode frequencies (slice [2..6])", pool.pool_id);
        print_top("  R1 (P7):", &freq_r1, 20);
        print_top("  R2 (P5):", &freq_r2, 20);
    }

    stats
}

/// Build a rust-bio FASTQ Reader from gzipped file
fn make_reader(path: &Path) -> FastqReader<BufReader<MultiGzDecoder<File>>> {
    let file = File::open(path).unwrap_or_else(|e| panic!("Failed to open {}: {}", path.display(), e));
    // Pass the *decoder* directly; Reader::new will wrap it in a BufReader internally.
    let gz = MultiGzDecoder::new(file);
    FastqReader::new(gz)
}

/// Build a rust-bio FASTQ Writer with gzip compression
fn make_writer(path: &Path) -> FastqWriter<GzEncoder<BufWriter<File>>> {
    let file = File::create(path).unwrap_or_else(|e| panic!("Failed to create {}: {}", path.display(), e));
    let enc = GzEncoder::new(BufWriter::new(file), Compression::default());
    FastqWriter::new(enc)
}

/// Canonicalize FASTQ IDs: remove trailing /1 or /2 and anything after first whitespace
fn canonical_id(id: &str) -> String {
    let head = id.split_whitespace().next().unwrap_or(id);
    if let Some(stripped) = head.strip_suffix("/1") {
        stripped.to_string()
    } else if let Some(stripped) = head.strip_suffix("/2") {
        stripped.to_string()
    } else {
        head.to_string()
    }
}

/// Uppercase a 4-mer slice (ASCII) into a Vec<u8>
fn upper4(s: &[u8]) -> Vec<u8> {
    s.iter()
        .map(|&b| if (b'a'..=b'z').contains(&b) { (b - b'a') + b'A' } else { b })
        .collect()
}

/// Reverse-complement a DNA slice (ASCII A/C/G/T/N). Returns uppercase Vec<u8>.
fn revcomp4(s: &[u8]) -> Vec<u8> {
    let map = |b: u8| match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    };
    let mut rc = Vec::with_capacity(4);
    for &b in s.iter().rev() {
        rc.push(map(b));
    }
    rc
}

/// Increment frequency map
fn bump(map: &mut HashMap<Vec<u8>, u64>, k: Vec<u8>) {
    *map.entry(k).or_insert(0) += 1;
}

/// Print top-k observed 4-mers
fn print_top(title: &str, map: &HashMap<Vec<u8>, u64>, k: usize) {
    let mut v: Vec<(&Vec<u8>, &u64)> = map.iter().collect();
    v.sort_by(|a, b| b.1.cmp(a.1));
    println!("{title}");
    for (i, (key, val)) in v.into_iter().take(k).enumerate() {
        println!("  {:>2}. {} : {}", i + 1, String::from_utf8_lossy(key), val);
    }
}