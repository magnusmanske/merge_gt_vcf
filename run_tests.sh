#!/usr/bin/env bash
\rm ./target/release/merge_gt_vcf
RUSTFLAGS="-C target-cpu=native" cargo build --release

CMDA=$(time ./target/release/merge_gt_vcf --serial < test_data/test.manifest | md5)
CMDB=$(cat test_data/test.md5)
if [ "$CMDA" = "$CMDB" ]; then
  echo "Plain VCF output OK"
else
  echo "Plain VCF output FAILED"
fi

CMDA=$(time ./target/release/merge_gt_vcf --serial --bgzip < test_data/test.manifest | md5)
CMDB=$(cat test_data/test.gz.md5)
if [ "$CMDA" = "$CMDB" ]; then
  echo "Bgzipped VCF output OK"
else
  echo "Bgzipped VCF output FAILED"
fi
