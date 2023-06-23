# Data preparation guide

The following commands show the steps taken to prepare the early and late datasets for this SARS-CoV-2 analysis.

``` bash
# Download full open metadata.
curl -OL https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst

# Select only high-quality sequences based on Nextclade
# QC status and select columns of interest.
zstd -d -c metadata.tsv.zst \
    | tsv-filter -H --str-eq QC_overall_status:good \
    | tsv-select -H -f strain,genbank_accession_rev,date,region,country,originating_lab,submitting_lab,Nextstrain_clade,Nextclade_pango \
    | zstd -c > filtered_metadata.tsv.zst

# Remove full metadata.
rm -f metadata.tsv.zst

# Force include the reference strain for rooting.
echo "Wuhan-Hu-1/2019" > reference_strain.txt

# Subsample early samples.
augur filter --metadata filtered_metadata.tsv.zst --include reference_strain.txt --min-date 2020-01-01 --max-date 2022-01-01 --group-by region week --subsample-max-sequences 2000 --output-strains early_global_strains.txt

# Assign random priorities with slight preference for
# recombinant lineages.
zstd -c -d filtered_metadata.tsv.zst | tsv-select -H -f strain,Nextclade_pango | sed 1d | awk 'OFS="\t" { priority = rand(); if (substr($2, 1, 1) == "X") { print $1,sprintf("%.2f", priority + 0.25) } else { print $1,sprintf("%.2f", priority) }}' > priorities.tsv

# Subsample late samples, prioritizing recombinants.
augur filter --metadata filtered_metadata.tsv.zst --include reference_strain.txt --min-date 2022-01-01 --group-by region week --subsample-max-sequences 2000 --priority priorities.tsv --output-strains late_global_strains.txt

# Download aligned sequences.
curl -OL https://data.nextstrain.org/files/ncov/open/aligned.fasta.zst

# Extract metadata and sequences for early and late samples.
augur filter --metadata filtered_metadata.tsv.zst --sequences aligned.fasta.zst --exclude-all --include early_global_strains.txt --output-metadata early_global_metadata.tsv --output-sequences early_global_aligned.fasta

augur filter --metadata filtered_metadata.tsv.zst --sequences aligned.fasta.zst --exclude-all --include late_global_strains.txt --output-metadata late_global_metadata.tsv --output-sequences late_global_aligned.fasta
```
