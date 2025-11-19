use crate::cigar::{clip, Cigar, CigarOp};
use itertools::Itertools;
use rust_htslib::bam::{self, IndexedReader, Record};

pub type Region<'a> = (&'a str, i64, i64);

#[derive(Debug)]
pub struct Prof {
    pub alts: Vec<f64>,
    pub depth: f64,
}

pub fn get_profile(bam: &mut IndexedReader, region: Region) -> Result<Prof, String> {
    let prof_len = (region.2 - region.1) as usize;
    let mut covs = vec![0; prof_len];
    let mut alts = vec![0; prof_len];
    bam.fetch(region).map_err(|e| e.to_string())?;

    let mut coverage = 0;
    for rec in bam::Read::records(bam) {
        let rec = rec.map_err(|e| e.to_string())?;

        if rec.is_secondary() || rec.is_supplementary() || rec.mapq() < 50 {
            continue;
        }
        update_profs(rec, &mut covs, &mut alts, region);
        coverage += 1;

        if coverage >= 200 {
            return Err("High depth".to_string());
        }
    }

    let depth = get_mean(&covs);

    let alts = alts
        .iter()
        .map(|v| *v as f64 / depth.max(1.0))
        .collect_vec();

    Ok(Prof { alts, depth })
}

fn update_profs(rec: Record, covs: &mut [u32], alts: &mut [u32], region: Region) {
    assert_eq!(covs.len() as i64, region.2 - region.1);
    assert_eq!(covs.len(), alts.len());
    let cigar = Cigar {
        ref_pos: rec.pos(),
        ops: rec.cigar().to_vec(),
    };
    let cigar = clip(&cigar, (region.1, region.2));
    if cigar.is_none() {
        return;
    }
    let cigar = cigar.unwrap();

    let mut ref_pos = cigar.ref_pos as usize;
    for op in cigar.ops {
        if ref_pos == region.2 as usize {
            break;
        }
        let index = ref_pos - region.1 as usize;
        match op {
            CigarOp::Match(len) | CigarOp::Equal(len) => {
                for cov in covs.iter_mut().skip(index).take(len as usize) {
                    *cov += 1;
                }
                ref_pos += len as usize;
            }
            CigarOp::Diff(len) => {
                for pos in index..index + len as usize {
                    alts[pos] += 1;
                    covs[pos] += 1;
                }
                ref_pos += len as usize;
            }
            CigarOp::Ins(len) => {
                alts[index] += len;
            }
            CigarOp::Del(len) => {
                for pos in index..index + len as usize {
                    alts[pos] += 1;
                    covs[pos] += 1;
                }
                ref_pos += len as usize;
            }
            CigarOp::SoftClip(_len) => {
                alts[index] += 1; // Avoid signal spikes from long and spurious softclips
            }
            CigarOp::HardClip(_) | CigarOp::Pad(_) | CigarOp::RefSkip(_) => {
                panic!("Missing logic to handle {:?}", op);
            }
        }
    }
}

fn get_mean(vals: &[u32]) -> f64 {
    vals.iter().sum::<u32>() as f64 / vals.len() as f64
}
