# Decoil

Decoil (deconvolve extrachromosomal circular DNA isoforms from long-read data) is a software package for reconstruction
circular DNA.

- [Getting started using docker or singularity](#gettingstarted)
- [Run example using docker or singularity](#testexample)
- [Run Decoil reconstruction using docker or singularity](#decoil-slim)
- [Decoil configurations](#decoil-config)
- [File formats](#decoil-file)
- [Citation](#citation)
- [License](#license)

## Getting started using docker or singularity<a name="gettingstarted"></a> 

As a prequisite you need to have install `docker` or `singularity` (you can install this from the official website or using `conda`).

### 1.1 Download as docker image

Download `decoil` docker image from `docker-hub`. This contains all the dependencies needed to run the software.<br/>
No additional installation needed. All the environment, packages, dependencies are all specified in the docker/singularity image. 


```commandline
# docker
docker pull madagiurgiu25/decoil:1.1.2-slim
```

### 1.2 Download as singularity image

```commandline
# singularity
singularity pull decoil.sif  docker://madagiurgiu25/decoil:1.1.2-slim
```

<br/>

### 2. Run example using docker or singularity <a name="testexample"></a> 

To test your installation check the [Example](docs/example.md).

<br/>

### 3. Run Decoil reconstruction using docker or singularity <a name="decoil-slim"></a> 

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
# docker
docker run -it --platform=linux/amd64 \
    -v ${BAM_INPUT}:/data/input.bam \
    -v ${BAM_INPUT}.bai:/data/input.bam.bai \
    -v ${GENOME}:/annotation/reference.fa \
    -v ${ANNO}:/annotation/anno.gtf \
    -v ${OUTPUT_FOLDER}:/mnt \
    -t madagiurgiu25/decoil:1.1.2-slim \
    decoil-pipeline sv-reconstruct \
            -b /data/input.bam \
            -r /annotation/reference.fa \
            -g /annotation/anno.gtf \
            -o /mnt -n ${NAME}
```

```bash
# singularity
mkdir -p ${OUTPUT_FOLDER}
singularity run \
    --bind ${BAM_INPUT}:/data/input.bam \
    --bind ${BAM_INPUT}.bai:/data/input.bam.bai \
    --bind ${GENOME}:/annotation/reference.fa \
    --bind ${ANNO}:/annotation/anno.gtf \
    --bind ${OUTPUT_FOLDER}:/mnt \
    decoil.sif \
    decoil-pipeline sv-reconstruct \
            -b /data/input.bam \
            -r /annotation/reference.fa \
            -g /annotation/anno.gtf \
            -o /mnt -n ${NAME}
```

<br/>

## Decoil configurations <a name="decoil-config"></a> 

An overview about the available functionalities:

|                	| [decoil-pipeline](#decoil-pipeline) 	| [decoil](#decoil-docs)	| [decoil-viz](#decoil-viz) 	|
|----------------	|--------	|-----------------	|------------	|
| SV calling     	| x       |                	 |            	|
| reconstruction 	| x      	| x               	|            	|
| visualization  	|        	|                 	| x          	|
| docker         	| x      	| x               	| x          	|
| singularity    	| x      	| x               	| x          	|

<br/>

### <a name="decoil-pipeline"></a> 1. Reconstruct ecDNA using `decoil-pipeline`

To reconstruct ecDNA we recommend to use `decoil-pipeline` using the `sv-reconstruct` mode.<br/>
This requires only a `.bam` file as input and generates internally all the files required for the reconstruction.


```commandline
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

The pipeline has the following [running modes](docs/decoil_pipeline_modes.md):

- `sv-only`
- `sv-reconstruct`
- `reconstruct-only`

<br/>

### <a name="decoil-viz"></a> 2. Visualization of ecDNA threads using `decoil-viz`

To interpret and visualize the results of the ecDNA reconstruction threads, use [decoil-viz](https://github.com/madagiurgiu25/decoil-viz).

<br/>

### <a name="decoil-docs"></a> 3. Advanced users only: reconstruct ecDNA using `decoil`

This configuration is the most flexible and allows users to use their own SV calls. For details go [here](docs/decoil_reconstruct.md).

<br/>

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

<br/>

## Citation <a name="citation"></a>

If you use Decoil for your work please cite our pre-print:

Madalina Giurgiu, Nadine Wittstruck, Elias Rodriguez-Fos, Rocio Chamorro Gonzalez, Lotte Bruckner, Annabell Krienelke-Szymansky, Konstantin Helmsauer, Anne Hartebrodt, Richard P. Koche, Kerstin Haase, Knut Reinert, Anton G. Henssen.
_Decoil: Reconstructing extrachromosomal DNA structural heterogeneity from long-read sequencing data_. bioRxiv, 2023, DOI: [https://doi.org/10.1101/2023.11.15.567169](https://www.biorxiv.org/content/10.1101/2023.11.15.567169v1)

## License <a name="license"></a> 

Decoil is distributed under the BSD 3-Clause license.  Consult the accompanying [LICENSE](LICENSE) file for more details.

## Disclaimer

Decoil and the content of this research-repository (i) is not suitable for a medical device; and (ii) is not intended
for clinical use of any kind, including but not limited to diagnosis or prognosis.
