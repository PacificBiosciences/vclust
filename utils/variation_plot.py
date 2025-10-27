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

DECODE_OP = {
    0: "M",  # match
    1: "I",  # insertion
    2: "D",  # deletion
    3: "N",  # refskip
    4: "S",  # softclip
    5: "H",  # hardclip
    6: "P",  # padding
    7: "=",  # match
    8: "X",  # mismatch
}

# Operations M, D, N, =, X consume reference
CONSUMES_REF = {0, 2, 3, 7, 8}

Profile = namedtuple("Profile", "covs alts")


def get_read_profile(read, locus_start, locus_end):
    ref_pos = read.reference_start
    covs = np.zeros(locus_end - locus_start)
    alts = np.zeros(locus_end - locus_start)

    for op_code, op_len in read.cigartuples:
        op_char = DECODE_OP.get(op_code, f"?{op_code}")

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


def plot_profiles(region_id, profiles, markers, image_path):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))
    fig.suptitle(f"Region: {region_id}", fontsize=12)

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
    region = args.region[0]
    chrom, start, end = region.replace(":", "-").split("-")
    start, end = int(start), int(end)

    pad = args.padding[0]
    profiles = []
    for reads_path in args.reads:
        if not os.path.exists(reads_path):
            print(f"Error: File not found -> {reads_path}", file=sys.stderr)
            continue
        profile = get_sample_profile(reads_path, chrom, start - pad, end + pad)
        profiles.append(profile)

    alts = np.array([prof.alts for prof in profiles])
    markers = [pad, end - start + pad]
    plot_profiles(region, alts, markers, args.image[0])


if __name__ == "__main__":
    main()
