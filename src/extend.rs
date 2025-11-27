use crate::locus::Locus;
use crate::models::{MODEL_REF, MODEL_VC, RADIUS};
use crate::profile::{get_profile, Prof};
use itertools::Itertools;
use logaddexp::LogAddExp;
use rust_htslib::bam::IndexedReader;

pub fn get_extension_offsets(
    locus: &Locus,
    bams: &mut Vec<IndexedReader>,
) -> Result<(i64, i64), String> {
    let region = extend_region(locus)?;

    let mut profs = Vec::new();
    for bam in bams {
        let prof = get_profile(bam, region)?;
        profs.push(prof);
    }

    let prof = merge(&profs);

    if prof.depth < 5.0 || prof.depth > 150.0 {
        return Err("DEPTH".to_string());
    }

    let alts = discretize(&prof.alts);

    let span = (RADIUS, RADIUS + locus.end - locus.start);
    let span =
        extend_to_ref_flanks(&alts, span, 150, locus.prior_vc).ok_or("EXTENSION".to_string())?;
    let span =
        extend_to_ref_flanks(&alts, span, 50, locus.prior_vc).ok_or("EXTENSION".to_string())?;
    let span =
        extend_to_ref_flanks(&alts, span, 25, locus.prior_vc).ok_or("EXTENSION".to_string())?;
    let span =
        extend_to_ref_flanks(&alts, span, 10, locus.prior_vc).ok_or("EXTENSION".to_string())?;

    let lf_offset = RADIUS - span.0;
    let rf_offset = span.1 - (RADIUS + locus.end - locus.start);

    Ok((lf_offset, rf_offset))
}

fn extend_region(locus: &Locus) -> Result<(&str, i64, i64), String> {
    if locus.start < RADIUS {
        Err("EXTENSION".to_string())
    } else {
        Ok((&locus.chrom[..], locus.start - RADIUS, locus.end + RADIUS))
    }
}

fn merge(profs: &[Prof]) -> Prof {
    let mut alts = Vec::new();
    for pos in 0..profs.first().unwrap().alts.len() {
        let alt = profs.iter().map(|prof| prof.alts[pos]).sum::<f64>() / profs.len() as f64;
        alts.push(alt);
    }

    let depth = profs.iter().map(|prof| prof.depth).sum::<f64>() / profs.len() as f64;

    Prof { alts, depth }
}

fn discretize(vals: &[f64]) -> Vec<u8> {
    vals.iter()
        .map(|val| {
            if *val <= 0.10 {
                0
            } else if *val < 0.25 {
                1
            } else if *val < 0.75 {
                2
            } else if *val < 1.50 {
                3
            } else if *val < 5.00 {
                4
            } else {
                5
            }
        })
        .collect()
}

fn extend_to_ref_flanks(
    alts: &[u8],
    span: (i64, i64),
    window_len: i64,
    prior_vc: f64,
) -> Option<(i64, i64)> {
    let mut lf_pos = span.0 - window_len;
    while lf_pos >= 0 {
        let window = &alts[lf_pos as usize..(lf_pos + window_len) as usize];
        let window = window.iter().rev().copied().collect_vec();
        let prob_ref = assess_window(&window[..], prior_vc);
        if prob_ref >= 0.5 {
            break;
        }
        lf_pos -= 1;
    }

    if lf_pos == 0 {
        return None;
    }

    let mut rf_pos = span.1;
    while rf_pos <= alts.len() as i64 - window_len {
        let window = &alts[rf_pos as usize..(rf_pos + window_len) as usize];
        let prob_ref = assess_window(window, prior_vc);

        if prob_ref >= 0.5 {
            break;
        }
        rf_pos += 1;
    }

    if alts.len() as i64 - window_len < rf_pos {
        return None;
    }

    Some((lf_pos + window_len, rf_pos))
}

fn assess_window(vals: &[u8], prior_vc: f64) -> f64 {
    let ll_norm = get_loglik(vals, &MODEL_REF) + (1.0 - prior_vc).ln();
    let ll_poly = get_loglik(vals, &MODEL_VC) + prior_vc.ln();
    let ll_sum = ll_norm.ln_add_exp(ll_poly);

    (ll_norm - ll_sum).exp()
}

fn get_loglik(prof: &[u8], model: &[f64; 1500]) -> f64 {
    let mut ll = 0.0;
    for (pos, val) in prof.iter().enumerate() {
        ll += model[pos * 6 + *val as usize].ln();
    }
    ll
}
