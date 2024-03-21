# Decoil-pipeline

To reconstruct ecDNA we recommend to use `decoil-pipeline` using the `sv-reconstruct` mode.<br/>
This requires only a `.bam` file as input and generates internally all the files required for the reconstruction.

The pipeline has the following running modes:

- [`sv-reconstruct`](#sv-reconstruct)
- [`sv-only`](#sv-only)
- [`reconstruct-only`](#reconstruct-only)

<br/>

```bash
# call help
docker run -it --platform=linux/amd64 -t madagiurgiu25/decoil:1.1.2-slim decoil-pipeline --help

usage: decoil-pipeline <workflow> <parameters> [<target>]
Example: 
    # run decoil including the processing and visualization steps
    decoil-pipeline -f sv-recontruct --bam <input> --outputdir <outputdir> --name <sample> --sv-caller <sniffles> -r <reference-genome> -g <annotation-gtf>
        

Decoil 1.1.2: reconstruct ecDNA from long-read data

positional arguments:
  {sv-only,sv-reconstruct,reconstruct-only}
                        sub-command help
    sv-only             Perform preprocessing
    sv-reconstruct      Perform preprocessing and reconstruction

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -n, --dry-run
  -f, --force
  -c, --use-conda
```

## `decoil-pipeline sv-reconstruct` mode: <a name="sv-reconstruct"></a> 

- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage
- reconstruct ecDNA using Decoil

Help:

```bash
# call help
docker run -it --platform=linux/amd64 -t madagiurgiu25/decoil:1.1.2-slim decoil-pipeline sv-reconstruct --help
```

Usage:

```bash
usage: decoil-pipeline sv-reconstruct --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>

Perform preprocessing (sv calling + coverage track) and reconstruction

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)
  -g ANNOTATION_GTF, --annotation-gtf ANNOTATION_GTF
                        GTF annotation
  -d, --debug           Debug mode
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
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated proportions (default: 0)

required named arguments:
  -b BAM, --bam BAM
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)

```

## `decoil-pipeline sv-only` mode: <a name="sv-only"></a> 

- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage

Help:

```bash
# call help
docker run -it --platform=linux/amd64 -t madagiurgiu25/decoil:1.1.2-slim decoil-pipeline sv-only --help
```

Usage:

```bash
usage: decoil-pipeline sv-only --bam <input> --outputdir <outputdir> --name <sample>

Perform preprocessing (sv calling + coverage track)

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)

required named arguments:
  -b BAM, --bam BAM
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
```

## `decoil-pipeline reconstruct-only` mode: <a name="reconstruct-only"></a> 

- reconstruct ecDNA using Decoil

Help:

```bash
# call help
docker run -it --platform=linux/amd64 -t madagiurgiu25/decoil:1.1.2-slim decoil-pipeline reconstruct-only --help 
```

Usage:

```bash
usage: decoil-pipeline reconstruct-only --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>

Perform reconstruction only

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads (default: 4)
  -g ANNOTATION_GTF, --annotation-gtf ANNOTATION_GTF
                        GTF annotation
  -d, --debug           Debug mode
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
  --filter-score FILTER_SCORE
                        Filter circular structures by estimated proportions (default: 0)

required named arguments:
  -b BAM, --bam BAM     Bam file
  -o OUTPUTDIR, --outputdir OUTPUTDIR
  --name NAME           Name of the sample
  -r REFERENCE_GENOME, --reference-genome REFERENCE_GENOME
                        Reference genome (fasta)
```