use chrono::Datelike;
use clap::Parser;
use locus::load_loci;
use rust_htslib::bam::IndexedReader;
use std::path::{Path, PathBuf};
use workflow::run_workflow;

mod cigar;
mod extend;
mod locus;
mod models;
mod profile;
mod workflow;

#[derive(Parser)]
#[command(name="HIFI-VCLUST",
          about="HiFi Variation Cluster Analysis Tool", 
          long_about = None,
          after_help = format!("Copyright (C) 2004-{} Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
          help_template = "{name} {version}\n{author}{about-section}\n{usage-heading}\n    {usage}\n\n{all-args}{after-help}",
          )]
#[command(arg_required_else_help(true))]
pub struct CliParams {
    #[clap(required = true)]
    #[clap(long = "genome")]
    #[clap(help = "Path to reference genome FASTA")]
    #[clap(value_name = "FASTA")]
    #[arg(value_parser = check_file_exists)]
    pub genome_path: PathBuf,

    #[clap(required = true)]
    #[clap(long = "reads")]
    #[clap(help = "BAM file with aligned HiFi reads")]
    #[clap(value_name = "READS")]
    #[clap(value_delimiter = ' ')]
    #[clap(num_args = 1..)]
    #[arg(value_parser = check_file_exists)]
    pub reads_paths: Vec<PathBuf>,

    #[clap(required = true)]
    #[clap(long = "regions")]
    #[clap(help = "BED file with region coordinates")]
    #[clap(value_name = "REGIONS")]
    #[arg(value_parser = check_file_exists)]
    pub repeats_path: PathBuf,
}

fn main() -> Result<(), String> {
    let args = CliParams::parse();
    let loci = load_loci(args.repeats_path)?;
    // TODO: Validate parameters

    let mut bams = Vec::new();
    for path in args.reads_paths {
        let bam = IndexedReader::from_path(&path).map_err(|e| e.to_string())?;
        bams.push(bam);
    }

    for locus in loci {
        if let Err(message) = run_workflow(&mut bams, &locus) {
            log::warn!("{message}");
        }
    }

    Ok(())
}

fn check_file_exists(path: &str) -> Result<PathBuf, String> {
    let path = Path::new(path);
    if path.exists() {
        Ok(path.to_path_buf())
    } else {
        Err(format!("File does not exist: {}", path.display()))
    }
}
