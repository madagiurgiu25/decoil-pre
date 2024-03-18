
You can test if Decoil is correctly installed by running the following example. 

## Download docker image

To run the example, download the `decoil:1.1.2-slim` docker image from `docker hub`.

```
# docker
docker pull madagiurgiu25/decoil:1.1.2-slim
```

## Download annotation data

```bash
# download annotation files
REFGENOME=reference.fa
GTFANNO=anno.gtf
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > $REFGENOME
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | gunzip -c > $GTFANNO
```


## Run decoil in `sv-reconstruct` mode

The example is started in `sv-reconstruct` mode and will perform:
- SV calling using sniffles1
- bigWig track generation using bamCoverage
- Decoil reconstruction

```bash

# run decoil in `sv-reconstruct` mode
docker run -it --platform=linux/amd64 \
    -v $PWD/test3:/examples \
    -v $PWD/$GTFANNO:/annotation/anno.gtf \
    -v $PWD/$REFGENOME:/annotation/reference.fa \
    -t madagiurgiu25/decoil:1.1.2-slim \
decoil-pipeline -f sv-reconstruct \
    --bam /examples/ecdna1/map.bam \
    --reference-genome /annotation/reference.fa \
    --annotation-gtf /annotation/anno.gtf \
    --outputdir /examples \
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