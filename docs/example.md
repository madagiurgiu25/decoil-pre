
# Example

You can test if Decoil is correctly installed by running the following example. 

## Download docker image and convert to singularity

If you want to run the example using `docker`, download the `decoil:1.1.2-slim` docker image from `docker hub`.

```bash
# docker
docker pull madagiurgiu25/decoil:1.1.2-slim
```

If you want to run the example using `singularity`, run the command below to convert the docker image into a singularity image file (SIF):

```bash
# singularity
singularity pull decoil.sif madagiurgiu25/decoil:1.1.2-slim
```

## Download annotation data

To run the example you need a reference genome and the genes annotation.

```bash
# download annotation files
REFGENOME=reference.fa
GTFANNO=anno.gtf
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > $REFGENOME
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | gunzip -c > $GTFANNO
```


## Run `decoil-pipeline` in `sv-reconstruct` mode

The example reconstructs ecDNA using `decoil-pipeline` in `sv-reconstruct` mode. This will perform:
- SV calling using sniffles1
- bigWig track generation using bamCoverage
- Decoil reconstruction

Using docker:

```bash
# run decoil-pipeline using docker `sv-reconstruct` mode
docker run -it --platform=linux/amd64 \
    -v $PWD/test3:/mnt \
    -v $PWD/$GTFANNO:/annotation/anno.gtf \
    -v $PWD/$REFGENOME:/annotation/reference.fa \
    -t madagiurgiu25/decoil:1.1.2-slim \
decoil-pipeline -f sv-reconstruct \
    --bam /examples/ecdna1/map.bam \
    --reference-genome /annotation/reference.fa \
    --annotation-gtf /annotation/anno.gtf \
    --outputdir /mnt \
    --name ecdna1
```

Using singularity:

```bash

# singularity needs to create upfront the output directory which will be mounted into the container
mkdir -p $PWD/test3

# run decoil-pipeline using singularity in `sv-reconstruct` mode
singularity run \
    --bind $PWD/test3:/mnt \
    --bind $PWD/$GTFANNO:/annotation/anno.gtf \
    --bind $PWD/$REFGENOME:/annotation/reference.fa \
    decoil.sif \
decoil-pipeline -f sv-reconstruct \
    --bam /examples/ecdna1/map.bam \
    --reference-genome /annotation/reference.fa \
    --annotation-gtf /annotation/anno.gtf \
    --outputdir /mnt \
    --name ecdna1
```

## Output

The command will create an output folder under `$PWD/test3`, containing the following files:

```commandline
tree $PWD/test3
|-- coverage.bw
|-- sv.sniffles.vcf
|-- sv.sniffles.bedpe
|-- logs
...
|-- reconstruct.ecDNA.bed
|-- reconstruct.ecDNA.filtered.bed
|-- reconstruct.ecDNA.fasta
|-- reconstruct.ecDNA.filtered.fasta
|-- reconstruct.links.ecDNA.txt
|-- reconstruct.linka.filtered.txt
|-- summary.txt
```