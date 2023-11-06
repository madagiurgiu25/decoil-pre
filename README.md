# Decoil

Decoil (deconvolve extrachromosomal circular DNA isoforms from long-read data) is a software package for reconstruction
circular DNA.

- [Getting started using docker](#gettingstarted)
- [Run Decoil reconstruction using docker](#decoil-slim)
- [Decoil configuration](#decoil-config)
- [File formats](#decoil-file)
- [Citation](#citation)
- [License](#license)

## Getting started using docker <a name="gettingstarted"></a> 

As a prequisite you need to have install `docker` (you can install this from the official website or using `conda`).
The advantage is that you do not need to install any environment as everything is packed inside the container.

### Download the docker image

This image contains all the dependencies needed to run the software.
No additional installation needed.


```commandline
# docker
docker pull madagiurgiu25/decoil:1.1.1-slim-test
```

#### Run example

You can test if Decoil is correctly installed by running the following example. This will create a output folder under `$PWD/test3`.
The example is started in `sv-reconstruct` mode which will perform:
- SV calling using sniffles1
- bigWig track generation using bamCoverage
- Decoil reconstruction

```bash
# download annotation files
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > $REFGENOME
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | gunzip -c > $GTFANNO

REFGENOME=reference.fa
GTFANNO=anno.gtf
# run decoil in `sv-reconstruct` mode
docker run -it --platform=linux/amd64 \
    -v $PWD/test3:/output \
    -v $PWD/$GTFANNO:/annotation/anno.gtf \
    -v $PWD/$REFGENOME:/annotation/reference.fa \
    -t madagiurgiu25/decoil:1.1.1-slim-test \
decoil -f sv-reconstruct \
    -b /examples/ecdna1/map.bam \
    -r /annotation/reference.fa \
    -g /annotation/anno.gtf \
    -o /output -n ecdna1
```

Check example output results:

```commandline
tree $PWD/test3
|-- clean.vcf
|-- config.json
|-- coverage.bw
|-- fragments_clean_1.bed
|-- fragments_clean_2.bed
|-- fragments_clean_3.bed
|-- fragments_initial.bed
|-- graph.gpickle
|-- logs
|   |-- logs_cov
|   |-- logs_decoil
|   |-- logs_sniffles
|   `-- logs_survivor
|-- metrics_frag_cov_distribution.png
|-- metrics_frag_len_cov_correlation.png
|-- metrics_frag_len_distribution.png
|-- reconstruct.bed
|-- reconstruct.bed_debug
|-- reconstruct.ecDNA.bed
|-- reconstruct.json
|-- reconstruct.links.txt_debug
|-- sv.sniffles.bedpe
`-- sv.sniffles.vcf
```

## Run Decoil reconstruction using docker <a name="decoil-slim"></a> 

To run Decoil on your data you need to cofigure the following parameters:

```bash
# run decoil with your input with standard parameters
BAM_INPUT="<absolute path to your BAM file>"
OUTPUT_FOLDER="<absolute path to your output folder>"
NAME="<sample name>"
GENOME="<absolute path to your reference genome file>"
ANNO="<absolute path to your gtf annotation file>"
```

and then run the following command:

```bash
docker run -it --platform=linux/amd64 \
    -v ${BAM_INPUT}:/data/input.bam \
    -v ${BAM_INPUT}.bai:/data/input.bam.bai \
    -v ${GENOME}:/annotation/reference.fa \
    -v ${ANNO}:/annotation/anno.gtf \
    -v ${OUTPUT_FOLDER}:/output \
    -t madagiurgiu25/decoil:1.1.1-slim-test \
    decoil -f sv-reconstruct \
            -b /data/input.bam \
            -r /annotation/reference.fa \
            -g /annotation/anno.gtf \
            -o /output -n ${NAME}
```

## Decoil configuration <a name="decoil-config"></a> 

To check Decoil's run modes:

```commandline
# call help
docker run -it --platform=linux/amd64 -t madagiurgiu25/decoil:1.1.1-slim-test decoil --help
```
```commandline
usage: decoil <workflow> <parameters> [<target>]
Example: 
    # run decoil including the processing and visualization steps
    decoil -f sv-recontruct --bam <input> --outputdir <outputdir> --name <sample> --sv-caller <sniffles> -r <reference-genome> -g <annotation-gtf>
        

Decoil 1.1.1: reconstruct ecDNA from long-read data

positional arguments:
  {sv-only,sv-reconstruct}
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

`sv-only` mode:
-
- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage

`sv-reconstruct` mode:
-
- SV calling using sniffles1
- convert VCF to BEDPE using SURVIVOR
- bigWig track generation using bamCoverage
- Decoil reconstruction

## File formats <a name="file-format"></a> 

The relevant output files for the users are:

- `reconstruct.bed` - contains all genomic fragments composing all reconstructions
- `reconstruct.ecDNA.bed` - contains all genomic fragments composing all the ecDNA labeled reconstructions
- `summary.txt` - summarize all the circular reconstructions

```commandline
cat reconstruct.bed

#chr    start   end     circ_id fragment_id     strand  coverage        estimated_proportions
chr2    15585356        15633376        0       5       +       149     75
chr3    11150000        11160001        0       41      -       103     75
chr3    11049997        11060001        0       33      +       117     75
chr2    15585356        15633376        3       5       +       149     36
chr3    11150000        11160001        3       41      -       103     36
chr3    11049997        11060001        3       33      +       117     36
chr2    15585356        15633376        3       5       +       149     36
chr2    16521052        16628305        3       13      +       37      36
chr3    10981202        11028470        3       25      -       31      36
chr12   68807722        68970910        2       53      +       252     252

```

```commandline
cat summary.txt

circ_id chr_origin      size(MB)        label   topology_idx    topology_name   estimated_proportions
0       chr3,chr2       0.068025                4       multi_region_inter_chr  75
3       chr3,chr2       0.270566        ecDNA   5       simple_duplications     36
2       chr12           0.163188        ecDNA   0       simple_circle           252
```

## Citation <a name="citation"></a>

## License <a name="license"></a> 

Decoil is distributed under [LICENSE](LICENSE), academic use-only.

 Disclaimer

Decoil and the content of this research-repository (i) is not suitable for a medical device; and (ii) is not intended
for clinical use of any kind, including but not limited to diagnosis or prognosis.
