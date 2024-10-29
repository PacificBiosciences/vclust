pub type CigarOp = rust_htslib::bam::record::Cigar;

#[derive(Debug, PartialEq, Clone)]
pub struct Cigar {
    pub ref_pos: i64,
    pub ops: Vec<CigarOp>,
}

/// Clips an alignment to a given reference region
///
/// Outputs reference start, query start, and operations of the clipped alignment
pub fn clip(cigar: &Cigar, region: (i64, i64)) -> Option<Cigar> {
    let (region_start, region_end) = region;
    let (read_start, read_end) = (cigar.ref_pos, get_reference_end(cigar));

    if read_end <= region_start || region_end <= read_start {
        return None;
    }

    let mut ref_pos = cigar.ref_pos;
    let mut op_iter = cigar.ops.iter();
    let mut current_op = op_iter.next();

    let mut clipped_ops = Vec::new();

    // Skip operations outside of the target region
    while current_op.is_some() && ref_pos + get_ref_len(current_op.unwrap()) <= region_start {
        ref_pos += get_ref_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    let mut clipped_ref_start = ref_pos;

    // Split operation overlapping the left flank (if any)
    if ref_pos < region_start {
        // Take care of situations where a single operation spans the entire region
        let ref_outside_len = region_start - ref_pos;
        let op_ref_len = get_ref_len(current_op.unwrap());
        let clipped_op_ref_len = if ref_pos + op_ref_len <= region_end {
            op_ref_len - ref_outside_len
        } else {
            region_end - region_start
        } as u32;
        clipped_ops.push(match current_op.unwrap() {
            CigarOp::Match(_) => CigarOp::Match(clipped_op_ref_len),
            CigarOp::RefSkip(_) => CigarOp::RefSkip(clipped_op_ref_len),
            CigarOp::Del(_) => CigarOp::Del(clipped_op_ref_len),
            CigarOp::Equal(_) => CigarOp::Equal(clipped_op_ref_len),
            CigarOp::Diff(_) => CigarOp::Diff(clipped_op_ref_len),
            op => panic!("Unexpected operation {:?}", op),
        });

        clipped_ref_start += ref_outside_len;

        ref_pos += get_ref_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    // Copy operations contained within the region
    while current_op.is_some() && ref_pos + get_ref_len(current_op.unwrap()) <= region_end {
        clipped_ops.push(*current_op.unwrap());
        ref_pos += get_ref_len(current_op.unwrap());
        current_op = op_iter.next();
    }

    // Split operation overlapping the right flank (if any)
    if current_op.is_some() && ref_pos < region_end {
        let ref_inside_len = region_end - ref_pos;
        clipped_ops.push(match current_op.unwrap() {
            CigarOp::Match(_) => CigarOp::Match(ref_inside_len as u32),
            CigarOp::RefSkip(_) => CigarOp::RefSkip(ref_inside_len as u32),
            CigarOp::Del(_) => CigarOp::Del(ref_inside_len as u32),
            CigarOp::Equal(_) => CigarOp::Equal(ref_inside_len as u32),
            CigarOp::Diff(_) => CigarOp::Diff(ref_inside_len as u32),
            op => panic!("Unexpected operation {:?}", op),
        });
    }

    let cigar = Cigar {
        ref_pos: clipped_ref_start,
        ops: clipped_ops,
    };
    Some(cigar)
}

fn get_reference_end(cigar: &Cigar) -> i64 {
    let mut ref_len = cigar.ref_pos;
    for op in &cigar.ops {
        ref_len += get_ref_len(op);
    }

    ref_len
}

fn get_ref_len(op: &CigarOp) -> i64 {
    match op {
        CigarOp::Match(len)
        | CigarOp::RefSkip(len)
        | CigarOp::Del(len)
        | CigarOp::Equal(len)
        | CigarOp::Diff(len) => *len as i64,
        CigarOp::Ins(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Pad(_) => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn make_cigar(ref_pos: i64, encoding: &str) -> Cigar {
        let ops = CigarString::try_from(encoding).unwrap().to_vec();
        Cigar { ref_pos, ops }
    }

    #[test]
    fn if_no_overlap_then_none() {
        //                     /AAATC
        //           CGC--TCGTTACG
        // CCCCCCCCCCCGCGGTCATTACGCCCCCCCCCC
        // |---10---||-----13----||---10---|
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let expected = clip(&cigar, (0, 10));
        assert_eq!(expected, None);

        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");
        let expected = clip(&cigar, (23, 33));
        assert_eq!(expected, None);
    }

    #[test]
    fn if_alignment_contained_inside_region_then_original_read() {
        let cigar = make_cigar(10, "5S3=2D2=1X2=5I3=10S");
        assert_eq!(clip(&cigar, (9, 23)), Some(cigar));
    }

    #[test]
    fn if_alignment_overlaps_left_flank_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");

        let expected = make_cigar(10, "3=2D");
        assert_eq!(clip(&cigar, (0, 15)), Some(expected));
    }

    #[test]
    fn if_op_overlaps_flanks_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");

        let expected = make_cigar(12, "1=2D2=");
        assert_eq!(clip(&cigar, (12, 17)), Some(expected));
    }

    #[test]
    fn if_op_spans_entire_region_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");

        let expected = make_cigar(21, "1=");
        assert_eq!(clip(&cigar, (21, 22)), Some(expected));
    }

    #[test]
    fn if_alignment_starts_inside_region_then_clipped_read() {
        // CGC--TCGTTAAATCACG <- Read
        // |||  || ||     |||
        // CGCAATCATT-----ACG <- Reference
        let cigar = make_cigar(10, "3=2D2=1X2=5I3=");

        let expected = make_cigar(10, "3=2D2=");
        assert_eq!(clip(&cigar, (0, 17)), Some(expected));
    }
}
