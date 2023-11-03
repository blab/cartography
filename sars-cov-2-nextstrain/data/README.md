# Data preparation guide

The following commands show the steps taken to prepare data for both early and late SARS-CoV-2 analysis.
These commands filter the full open metadata to high-quality records and the columns we need for this project.
The final step uploads the filtered metadata to an AWS S3 bucket with a date-stamp for future reproducibility.

``` bash
# Download full open metadata.
curl -OL https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst

# Select only high-quality sequences based on Nextclade
# QC status and select columns of interest.
zstd -d -c metadata.tsv.zst \
    | tsv-filter -H --str-eq QC_overall_status:good \
    | tsv-filter -H --or --str-eq "clock_deviation:?" --le clock_deviation:20 \
    | tsv-filter -H --or --str-eq "clock_deviation:?" --ge clock_deviation:-20 \
    | tsv-select -H -f strain,genbank_accession_rev,date,region,country,division,originating_lab,submitting_lab,Nextstrain_clade,Nextclade_pango \
    | zstd -c > filtered_metadata.tsv.zst

# Download full open aligned sequences.
curl -fsSLO --proto '=https' https://data.nextstrain.org/files/ncov/open/aligned.fasta.zst

# Get the list of strains with metadata.
zstd -c -d filtered_metadata.tsv.zst | tsv-select -H -f strain | sed 1d | sort -k 1,1 | uniq > filtered_strains.txt

# Get the sequences that match the metadata strains.
seqkit grep -f filtered_strains.txt aligned.fasta.zst | zstd -c > filtered_aligned.fasta.zst

# Upload date-stamped filtered metadata to AWS S3 bucket.
aws s3 cp \
    filtered_metadata.tsv.zst \
    s3://nextstrain-data/files/workflows/cartography/filtered_metadata_$(date "+%Y-%m-%d").tsv.zst

# Upload date-stamped filtered alignments to AWS S3 bucket.
aws s3 cp \
    filtered_aligned.fasta.zst \
    s3://nextstrain-data/files/workflows/cartography/filtered_aligned_$(date "+%Y-%m-%d").fasta.zst

# Remove full metadata and alignments.
rm -f metadata.tsv.zst aligned.fasta.zst
```
