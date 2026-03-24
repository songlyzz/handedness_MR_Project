#!/bin/bash

# directory
work_dir="./hand_LDSC_analysis"
ref_dir="./software/ldsc/reference"
soft_dir="./software/ldsc"
THREADS=8
Size=580148 # gwas population
cd ${work_dir}
source activate ldsc
# pre ldsc data
python ${soft_dir}/munge_sumstats.py --sumstats ${work_dir}/XX_munge.txt --N 580148 --out ${work_dir}/df_ldsc --chunksize 500000 --merge-alleles ${ref_dir}/eur_w_ld_chr/w_hm3.snplist

echo $(date +%Y-%m-%d-%H)
echo $(date +%Y-%m-%d-%H)
echo "============== ldsc format data already prepared =============="
echo "============== ldsc format data already prepared =============="

# run
# heritability
python ${soft_dir}/ldsc.py \
--h2 df_ldsc.sumstats.gz \
--ref-ld-chr ${ref_dir}/eur_w_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ \
--out df_h2

# correlation
python ${soft_dir}/ldsc.py \
--rg df1_ldsc.sumstats.gz,df2_ldsc.sumstats.gz \
--ref-ld-chr ${ref_dir}/eur_w_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ \
--out overlap

echo "============== ldsc analysis =============="
echo "============== ldsc analysis =============="
