# Decoil

For advanced users only: you can run `decoil reconstruct` to reconstruct ecDNA.<br/>
This requires as input the SV calls (.vcf), the aligment file (.bam) and the coverage track (.bw), which needs to be computed upfront by the user.


Usage:

```bash
usage: decoil reconstruct -i <vcffile> -c <coveragefile> -b [bamfile] --outputdir <outputdir> --name <sample> -r <reference_genome>

options:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Bam file (not required if --fast enabled)
  -g ANNOTATION_GTF, --annotation-gtf ANNOTATION_GTF
                        GTF annotation
  -d, --debug           Debug mode
  --extend-allowed-chr EXTEND_ALLOWED_CHR
                        Add list custom assemblies/chromosomes (e.g. "chr1,chr2,chr3")
  --fast                Reconstruct fast (not accurate and does not require a bam file)
  --multi               Multi-sample VCF file
  --skip                Skip creation of fasta files (to save space)
  --min-sv-len MIN_SV_LEN
                        Minimal SV length (default: 500X)
  --fragment-min-cov FRAGMENT_MIN_COV
                        Minimal fragment coverage (default: 5X)
  --fragment-max-cov FRAGMENT_MAX_COV
                        Maximal fragment coverage (default: 100000X)
  --fragment-min-size FRAGMENT_MIN_SIZE
                        Minimal fragment size (default: 500bp)
  --min-vaf MIN_VAF     Minimal VAF acceptance SV (default: 0.01)
  --min-cov-alt MIN_COV_ALT
                        Minimal supporting reads SV (default: 6X)
  --max-explog-threshold MAX_EXPLOG_THRESHOLD
                        Maximal score; better not change this (default: 0.1)
  --min-cov MIN_COV     Minimal coverage on site (default: 8X)
  --sv-caller SV_CALLER
                        SV caller name matching the VCF input ['sniffles1', 'sniffles2', 'cutesv', 'nanomonsv', 'lumpy', 'delly', 'multivcf']
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated copy-number (default: 0)

required named arguments:
  -c COVERAGE, --coverage COVERAGE
                        Coverage file (bigwig)
  -i VCF, --vcf VCF     Vcf file
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)
```
