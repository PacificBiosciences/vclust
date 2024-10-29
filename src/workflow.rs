use crate::extend::get_extension_offsets;
use crate::locus::Locus;
use rust_htslib::bam::IndexedReader;

pub fn run_workflow(bams: &mut Vec<IndexedReader>, locus: &Locus) -> Result<(), String> {
    let in_region = format!("{}:{}-{}", locus.chrom, locus.start, locus.end);
    let offsets = get_extension_offsets(locus, bams);

    if let Some((lf, rf)) = offsets {
        let out_region = format!("{}:{}-{}", locus.chrom, locus.start - lf, locus.end + rf);
        println!("{}\t{in_region}\t{lf}\t{rf}\t{out_region}", locus.name);
    } else {
        println!("{}\t{in_region}\tNA", locus.name);
    }

    Ok(())
}
