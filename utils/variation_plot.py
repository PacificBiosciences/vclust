# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "pysam",
# ]
# ///

import argparse
import os
import sys
import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from matplotlib.colors import LinearSegmentedColormap

# SAM/BAM CIGAR op codes -> characters (per hts-spec)
CIGAR_CODE_TO_CHAR = {
    0: "M",  # alignment match (can be = or X)
    1: "I",  # insertion to the reference
    2: "D",  # deletion from the reference
    3: "N",  # skipped region from the reference (e.g., intron)
    4: "S",  # soft clipping (clipped sequence present in SEQ)
    5: "H",  # hard clipping (clipped sequence NOT present in SEQ)
    6: "P",  # padding (silent deletion from padded reference)
    7: "=",  # sequence match
    8: "X",  # sequence mismatch
    # 9: "B"  # back (deprecated in SAM spec; not produced by modern tools)
}

# Which ops consume reference positions (advance the reference coordinate)
CONSUMES_REF = {0, 2, 3, 7, 8}  # M, D, N, =, X

Profile = namedtuple("Profile", "covs alts")


def get_read_profile(read, locus_start, locus_end):
    # print(locus_start, locus_end)
    # rname = bam.get_reference_name(read.reference_id)
    ref_pos = read.reference_start  # 0-based leftmost coordinate on reference
    covs = np.zeros(locus_end - locus_start)
    alts = np.zeros(locus_end - locus_start)

    # print(
    #    f"\nREAD: {read.query_name}  |  REF: {rname}  |  start={read.reference_start}"
    # )
    for op_code, op_len in read.cigartuples:
        op_char = CIGAR_CODE_TO_CHAR.get(op_code, f"?{op_code}")
        # print(op_char, op_len)

        if op_char in ["=", "M"]:
            op_start, op_end = ref_pos, ref_pos + op_len
            for op_pos in range(max(op_start, locus_start), min(op_end, locus_end)):
                assert op_pos - locus_start >= 0 and locus_end - op_pos >= 0
                covs[op_pos - locus_start] += 1
            ref_pos += op_len
        elif op_char == "I":
            if locus_start <= ref_pos < locus_end:
                alts[ref_pos - locus_start] += op_len
        elif op_char in ["D", "X"]:
            op_start, op_end = ref_pos, ref_pos + op_len
            for op_pos in range(max(op_start, locus_start), min(op_end, locus_end)):
                assert op_pos - locus_start >= 0 and locus_end - op_pos >= 0
                covs[op_pos - locus_start] += 1
                alts[op_pos - locus_start] += 1
            ref_pos += op_len
        elif op_char == "S":
            if locus_start <= ref_pos < locus_end:
                alts[ref_pos - locus_start] += 1
        else:
            assert False, f"Unknown op char {op_char}"

    return Profile(covs, alts)


def get_sample_profile(bam_path, chrom, start, end):
    """
    Load reads from a BAM/CRAM using pysam, parse each read's CIGAR, and print:
      - CIGAR operation character (e.g., M, I, D, N, S, H, =, X)
      - operation length
      - reference span [start, end) for ops that consume reference bases

    Parameters
    ----------
    bam_path : str
        Path to an indexed BAM/CRAM file (e.g., ".bam" with ".bai" index present).
    region : str | None
        Optional region string like "chr1:100000-101000". If None, iterates all reads.

    Notes
    -----
    - Reference coordinates are 0-based half-open [start, end).
    - For ops that do not consume reference (I, S, H, P), the reference position does not advance
      and we print the current reference coordinate (no span).
    - Unmapped reads or reads without a CIGAR are skipped.
    """
    sample_covs, sample_alts = np.zeros(end - start), np.zeros(end - start)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        iterator = bam.fetch(chrom, start, end)
        for read in iterator:
            if read.is_unmapped or read.cigartuples is None:
                continue
            if read.is_supplementary or read.is_secondary:
                continue
            if read.mapq < 50:
                continue

            read_profile = get_read_profile(read, start, end)
            sample_covs += read_profile.covs
            sample_alts += read_profile.alts
    mean_depth = max(1.0, np.mean(sample_covs))
    sample_alts /= mean_depth
    return Profile(sample_covs, sample_alts)


def plot_profiles(profiles, markers, image_path):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))
    locus = "Locus"
    fig.suptitle(f"{locus} locus", fontsize=12)

    colors = [(0.0, "#DDDDDD"), (1.0, "#7618DC")]
    cmap = LinearSegmentedColormap.from_list("gray_to_purple", colors, N=5)
    plot1 = ax.matshow(np.minimum(profiles, 1.0), aspect="auto", cmap=cmap)
    fig.colorbar(
        plot1,
        ax=ax,
        location="bottom",
        fraction=0.05,
        shrink=0.5,
        label="Fraction of alt bases",
    )
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Sample Index")
    ax.xaxis.set_ticks_position("bottom")

    for marker in markers:
        ax.axvline(x=marker, color="#e63946", linestyle="-")

    plt.tight_layout()
    plt.savefig(image_path)


def parse_args():
    parser = argparse.ArgumentParser(description="A tool to generate variation plots")
    parser.add_argument(
        "--reads", required=True, nargs="+", help="One or more BAM files"
    )
    parser.add_argument("--region", required=True, nargs=1, help="Region to analyze")
    parser.add_argument(
        "--padding",
        required=True,
        type=int,
        nargs=1,
        help="Padding around the region",
    )
    parser.add_argument(
        "--image",
        required=True,
        nargs=1,
        help="Output image path (must end in .pdf or .png)",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    pad = args.padding[0]
    profiles = []
    for reads_path in args.reads:
        if not os.path.exists(reads_path):
            print(f"Error: File not found -> {reads_path}", file=sys.stderr)
            continue
        chrom, start, end = args.region[0].replace(":", "-").split("-")
        start, end = int(start), int(end)
        profile = get_sample_profile(reads_path, chrom, start - pad, end + pad)
        profiles.append(profile)

    alts = np.array([prof.alts for prof in profiles])
    markers = [pad, end - start + pad]
    plot_profiles(alts, markers, args.image[0])

    # bam_files.append(pysam.AlignmentFile(reads_path))

    # bam_file = bam_files[0]
    # print(bam_file)

    # print(f"Loading file: {reads_path}")
    # with open(reads_path, "r") as f:
    #    # Example: print first 100 characters
    #    content = f.read(100)
    #    print(f"First 100 characters of {reads_path}:\n{content}\n---")


if __name__ == "__main__":
    main()
