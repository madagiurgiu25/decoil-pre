# Decoil

For advanced users only: you can run `decoil reconstruct` to reconstruct ecDNA.<br/>
This requires as input the SV calls (.vcf), the aligment file (.bam) and the coverage track (.bw), which needs to be computed upfront by the user.


Usage:

```bash
usage: decoil reconstruct -b <bamfile> -i <vcffile> -c <coveragefile> --outputdir <outputdir> --name <sample> -r <reference_genome>

optional arguments:
  -h, --help            show this help message and exit
  -g ANNOTATION_GTF, --annotation-gtf ANNOTATION_GTF
                        GTF annotation
  -d, --debug           Debug mode
  --fast                Reconstruct fast (not accurate and does not require a bam file)
  --min-sv-len MIN_SV_LEN
                        Minimal SV length (default: 500X)
  --fragment-min-cov FRAGMENT_MIN_COV
                        Minimal fragment coverage (default: 5X)
  --fragment-min-size FRAGMENT_MIN_SIZE
                        Minimal fragment size (default: 500bp)
  --min-vaf MIN_VAF     Minimal VAF acceptance SV (default: 0.01)
  --min-cov-alt MIN_COV_ALT
                        Minimal supporting reads SV (default: 6X)
  --max-explog-threshold MAX_EXPLOG_THRESHOLD
                        Maximal score; better not change this (default: 0.1)
  --min-cov MIN_COV     Minimal coverage on site (default: 8X)
  --sv-caller SV_CALLER
                        SV caller name {sniffles, sniffles2, cutesv}
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated copy-number (default: 0)

required named arguments:
  -b BAM, --bam BAM     Bam file
  -c COVERAGE, --coverage COVERAGE
                        Coverage file (bigwig)
  -i VCF, --vcf VCF     Vcf file
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)
```
