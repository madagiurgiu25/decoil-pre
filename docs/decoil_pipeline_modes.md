# Decoil-pipeline

To reconstruct ecDNA we recommend to use `decoil-pipeline` using the `sv-reconstruct` mode.<br/>
This requires only a `.bam` file as input and generates internally all the files required for the reconstruction.

The pipeline has the following running modes:

- [`sv-only`](#sv-only)
- [`sv-reconstruct`](#sv-reconstruct)
- [`reconstruct-only`](#reconstruct-only)

<br/>

```bash
# call help
usage:
                  decoil-pipeline [options] <run-mode> <parameters> [<target>]
                  or
                  decoil <run-mode> <paramters>

Decoil 2.0.0a1: reconstruct ecDNA from long-read data

positional arguments:
  {sv-only,sv-reconstruct,reconstruct-only}
                        sub-command help

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -n, --dry-run
  -f, --force
```

## `decoil-pipeline sv-only` mode: <a name="sv-only"></a> 

- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage

Help:

```bash
# call help
usage: decoil-pipeline sv-only --bam <input> --outputdir <outputdir> --name <sample>

Perform preprocessing (sv calling + coverage track)

options:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)
  --sv-caller SV_CALLER
                        SV caller name matching the VCF input ['sniffles1', 'sniffles2']

required named arguments:
  -b BAM, --bam BAM
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
```


## `decoil-pipeline sv-reconstruct` mode: <a name="sv-reconstruct"></a> 

- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage
- reconstruct ecDNA using Decoil

Help:

```bash
# call help
usage: decoil-pipeline sv-reconstruct --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>

Perform preprocessing (sv calling + coverage track) and reconstruction

options:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)
  -g ANNOTATION_GTF, --annotation-gtf ANNOTATION_GTF
                        GTF annotation
  -d, --debug           Debug mode
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
                        SV caller name matching the VCF input ['sniffles1', 'sniffles2']
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated proportions (default: 0)
  --extend-allowed-chr EXTEND_ALLOWED_CHR
                        Add list custom assemblies/chromosomes (e.g. --extend-allowed-chr chr1,chr2,chr3
  --skip                Skip creation of fasta files (to save space)
  --filt-version FILT_VERSION
                        Sniffles2 filter version

required named arguments:
  -b BAM, --bam BAM
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)
```

## `decoil-pipeline reconstruct-only` mode: <a name="reconstruct-only"></a> 

- reconstruct ecDNA using Decoil

Help:

```bash
# call help
usage: decoil-pipeline reconstruct-only --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>

Perform reconstruction only

options:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)
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
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated proportions (default: 0)
  --sv-caller SV_CALLER
                        SV caller name matching the VCF input ['sniffles1', 'sniffles2', 'cutesv', 'nanomonsv', 'lumpy', 'delly', 'multivcf']

required named arguments:
  -b BAM, --bam BAM     Bam file
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)
```