<h1 align="center">vclust</h1>

<p align="center">A tool to detect variation clusters in PacBio HiFi sequencing data</p>

***

Variation clusters are loosely defined as contiguous regions that contain some
variation in a given set of genomes. The `vclust` method is designed to locate
variation clusters around an input set of seed regions using variation extracted
from PacBio HiFi whole-genome sequencing data.

Authors: [Egor Dolzhenko](mailto:edolzhenko@pacificbiosciences.com)

## Installation and usage

`vclust` can be built from source (requires Rust, which can be installed with the
[rustup](https://rustup.rs/) tool).

```bash
git clone git@github.com:PacificBiosciences/vclust.git
cd vclust
cargo build --release
```

Once installed, `vclust` can be run like so:

```bash
./vclust --genome genome.fa --reads bams/*.bam --regions regions.bed > extended_regions.txt
```

where

* `genome.fa` is the reference genome (the same reference genome as used for
    read alignment)
* `bams/*.bam` paths to the aligned PacBio HiFi BAM files
* `regions.bed` a BED file with seed regions to extend
* `extended_regions.txt` an output file with the extended regions

The input file `regions.bed` is expected to contain coordinates and identifiers
of the regions to be profiled.

```csv
chr1    57367043        57367119        region1
chr1    146228800       146228821       region2
chr1    149390802       149390841       region3
chr10   79826383        79826404        region4
chr10   93702522        93702547        region5
chr11   66744821        66744850        region6
chr11   119206289       119206322       region7
```

The output file `extended_regions.txt` contains region identifiers (column 1),
coordinates of the original input regions (column 2), input regions' start / end
extension lengths (columns 3,4), and coordinates of the corresponding extended
regions (column 5). If `vclust` is unable to extend a given region (due to, say,
the lack of read coverage), the third column is set to `NA`.

```csv
region1 chr1:57367043-57367119      0    0       chr1:57367043-57367119
region2 chr1:146228800-146228821   NA
region3 chr1:149390802-149390841    0    0       chr1:149390802-149390841
region4 chr10:79826383-79826404    72    0       chr10:79826311-79826404
region5 chr10:93702522-93702547     0    0       chr10:93702522-93702547
region6 chr11:66744821-66744850     0    0       chr11:66744821-66744850
region7 chr11:119206289-119206322   0    0       chr11:119206289-119206322
```

In the example above, `vclust` did not identify significant variation around
`region1` and hence the extension lengths were set to 0. The reported extended
region in column 5 is the original region itself. On the other hand, the start
of `region4` was extended by 72 bps, expanding coordinates of the original
region from `chr10:79826383-79826404` to `chr10:79826311-79826404`.

## Citation

vclust is described in the Methods section of this paper:

[Ben Weisburd, Egor Dolzhenko,  Mark F. Bennett, Matt C. Danzi, Adam English,  Laurel Hiatt, Hope Tanudisastro, Nehir Edibe Kurtas, Helyaneh Ziaei Jam, Harrison Brand, Fritz J. Sedlazeck, Melissa Gymrek, Harriet Dashnow,  Michael A. Eberle, Heidi L. Rehm. Defining a tandem repeat catalog and variation clusters for genome-wide analyses and population databases. _bioRxiv_, 2024-10.](https://www.biorxiv.org/content/10.1101/2024.10.04.615514v1)

## Need help?

If you notice any missing features, bugs, or need assistance with analyzing the output
of vclust, please don't hesitate to open a GitHub issue or contact the developers by
[email](mailto:edolzhenko@pacificbiosciences.com).

vclust is a pre-release software intended for research use only and not for use in
diagnostic procedures. While efforts have been made to ensure that vclust lives up
to the quality that PacBio strives for, we make no warranty regarding this software.

As vclust is not covered by any service level agreement or the like, please do not
contact a PacBio Field Applications Scientists or PacBio Customer Service for
assistance with any vclust release. Please report all issues through GitHub instead.
We make no warranty that any such issue will be addressed, to any extent or within
any time frame.

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
