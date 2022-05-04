# merge_gt_vcf
A tool to merge VCF files that have the exact same number of (data) rows, chromosomes, positions, and ref/alt alleles. It takes a list of "genotying" VCF files from STDIN and outputs plain-text VCF, consisting of the first VCF file in the list, with sample columns from all the other VCF files appended.

**Caution**: This _only_ works if all VCF files have the exact same `CHROM`, `POS`, `ID`, `REF`, and `ALT` columns. There is a check for this in the code but currently commented out for performance reasons. It _will_ fail if the number of data rows is not exactly the same in all files.

# Installation
1. [Install Rust](https://www.rust-lang.org/tools/install).
2. `git clone https://github.com/magnusmanske/merge_gt_vcf.git`
3. `cargo build --release`

# Usage
```
# For plain-text VCF:
./merge_gt_vcf < FILE_WITH_VCF_FILENAMES > OUTPUT.vcf

# For bgzipped VCF:
./merge_gt_vcf < FILE_WITH_VCF_FILENAMES | bgzip > OUTPUT.vcf.gz
```
