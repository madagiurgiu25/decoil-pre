# Decoil

Decoil (deconvolve extrachromosomal circular DNA isoforms from long-read data) is a software package for reconstruction
circular DNA.

- [Getting started using docker](#gettingstarted)
- [Citation](#citation)
- [License](#license)

### Getting started using docker <a name="gettingstarted"></a> 

As a prequisite you need to have install `docker` (you can install this from the official website or using `conda`).

#### Download the image

```commandline
# docker
docker pull decoil:1.1.0beta5
```

#### Run example

You can test if Decoil is correctly installed by running the following example.

```bash
# download annotation files
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz | gunzip -c > $REFGENOME
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.basic.annotation.gtf.gz | gunzip -c > $GTFANNO

REFGENOME=reference.fa
GTFANNO=anno.gtf

# run decoil in `reconstruct` mode
docker run -v $PWD/testdecoil:/output \
-v $PWD/$GTFANNO:/annotation/anno.gtf \
-v $PWD/$REFGENOME:/annotation/reference.fa \
-t decoil:1.1.0beta5 \
decoil reconstruct \
-i /decoil/examples/ecdna1/sv.sniffles.vcf \
-b /decoil/examples/ecdna1/coverage.bw \
-x /decoil/examples/ecdna1/map.bam \
-r /annotation/reference.fa \
-g /annotation/anno.gtf \
-o /output -n ecdna1 --plot True
```

### File formats

### Citation <a name="citation"></a>

### License <a name="license"></a> 

Decoil is distributed under [LICENSE](LICENSE), academic use-only.

### Disclaimer

Decoil and the content of this research-repository (i) is not suitable for a medical device; and (ii) is not intended
for clinical use of any kind, including but not limited to diagnosis or prognosis.
