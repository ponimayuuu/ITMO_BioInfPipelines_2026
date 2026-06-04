# Instrucitons

1. To run this pipeline propely firdtly install nf-core bcftools module in a pipline directory

 ```nf-core modules install bcftools/mpileup```

2. Download samples

```for SRR in SRR13191702 SRR13191703 SRR14031324 SRR14031325; do fasterq-dump --split-files $SRR gzip ${SRR}_1.fastq ${SRR}_2.fastq done```

3. Download reference

```wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -O sars_cov2_ref.fasta```     

4. Enjoy!

 ```nextflow run main.nf -profile local --sample_csv samples.csv --reference sars_cov2_ref.fasta```
