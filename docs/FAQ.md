# FAQ

### I have a CRAM file instead of a BAM file. What should I do?

Convert the CRAM file to BAM before running Decoil. In most cases, working with BAM will be faster and simpler.

```
samtools view -T reference.fa -b -o output.bam input.cram
samtools index output.bam
```


---

### I get a `remove_duplicated_edges` error. What should I do?

This issue should be fixed in newer versions of Decoil. Please upgrade to the latest stable release:

```bash
python -m pip install decoil==2.0.0
```

---

### My sample runs for a very long time. What should I do?

Some samples contain regions with extremely high coverage. These are not necessarily focal amplifications and may instead result from mapping artifacts or repetitive regions. Very high local coverage can also confuse the SV caller and substantially increase runtime.

You can try the following:

1. Ignore genomic bins with extremely high coverage:

```bash
--fragment-max-cov 10000
```

This excludes bins with coverage above 10,000×.

2. Downsample the BAM file to approximately **5× mean coverage**.

3. Use more stringent SV filtering. For example, for WGS data sequenced at approximately **30× mean coverage**, you can try:

```bash
--min-vaf 0.1
--min-cov-alt 30
--fragment-min-cov 50
--min-cov 30
--min-sv-len 1000
```

These values are intended as a practical starting point and may need adjustment depending on the coverage distribution and quality of your sample.

---

### How can I use Decoil with short-read sequencing?

Decoil supports a fast-mode that enables reconstructions from short-reads. For this mode only SV.vcf and Coverage.bw are required, no BAM files needed.

```
decoil reconstruct <other params> --fast --sv-caller lumpy
```

### Can I use Decoil with a 

### Can I use Decoil with a custom genome?

Yes you can. Please provide (1) a `.fasta` while running decoil and (2) include all the non-canonical chromosomes sequences names. This generalize to viral or bacterial genomes. 

```
decoil reconstruct <other params> -r custom.fasta --extend-allowed-chr "chrCustom1,chrCustom2"
```

If you intend to visualize your structure with [`decoil-viz`](https://github.com/madagiurgiu25/decoil-viz) you will also need matched `.gtf`.

### Does Decoil support T2T?

Yes it does. You need to provide the T2T `.fasta` as parameter `-r T2T.fasta. You do not need to specify the chromosome names using `--allowed-chr` as the names of the canonical chromosomes in the UCSC, NCBI and EMBL annotation for hg19, hg38 and T2T are included.

See below the list of chromosomal names included by default when using Decoil:

```
ALLOWED_CHR = [
		"1",
		"2",
		"3",
		"4",
		"5",
		"6",
		"7",
		"8",
		"9",
		"10",
		"11",
		"12",
		"13",
		"14",
		"15",
		"16",
		"17",
		"18",
		"19",
		"20",
		"21",
		"22",
		"X",
		"Y",
		"M",
		"chr1",
		"chr2",
		"chr3",
		"chr4",
		"chr5",
		"chr6",
		"chr7",
		"chr8",
		"chr9",
		"chr10",
		"chr11",
		"chr12",
		"chr13",
		"chr14",
		"chr15",
		"chr16",
		"chr17",
		"chr18",
		"chr19",
		"chr20",
		"chr21",
		"chr22",
		"chrX",
		"chrY",
		"chrM",
		"CP068255.2",
		"CP068256.2",
		"CP068257.2",
		"CP068258.2",
		"CP068259.2",
		"CP068260.2",
		"CP068261.2",
		"CP068262.2",
		"CP068263.2",
		"CP068264.2",
		"CP068265.2",
		"CP068266.2",
		"CP068267.2",
		"CP068268.2",
		"CP068269.2",
		"CP068270.2",
		"CP068271.2",
		"CP068272.2",
		"CP068273.2",
		"CP068274.2",
		"CP068275.2",
		"CP068276.2",
		"CP068277.2",
		"CP086569.2",
		"CP068254.1",
		"NC_060925.1",
		"NC_060926.1",
		"NC_060927.1",
		"NC_060928.1",
		"NC_060929.1",
		"NC_060930.1",
		"NC_060931.1",
		"NC_060932.1",
		"NC_060933.1",
		"NC_060934.1",
		"NC_060935.1",
		"NC_060936.1",
		"NC_060937.1",
		"NC_060938.1",
		"NC_060939.1",
		"NC_060940.1",
		"NC_060941.1",
		"NC_060942.1",
		"NC_060943.1",
		"NC_060944.1",
		"NC_060945.1",
		"NC_060946.1",
		"NC_060947.1",
		"NC_060948.1",
		"NC_000001.11",
	"NC_000002.12",
	"NC_000003.12",
	"NC_000004.12",
	"NC_000005.10",
	"NC_000006.12",
	"NC_000007.14",
	"NC_000008.11",
	"NC_000009.12",
	"NC_000010.11",
	"NC_000011.10",
	"NC_000012.12",
	"NC_000013.11",
	"NC_000014.9",
	"NC_000015.10",
	"NC_000016.10",
	"NC_000017.11",
	"NC_000018.10",
	"NC_000019.10",
	"NC_000020.11",
	"NC_000021.9",
	"NC_000022.11",
	"NC_000023.11",
	"NC_000024.10",
	"NC_012920.1",
	"NC_000067.7",
	"NC_000068.8",
	"NC_000069.7",
	"NC_000070.7",
	"NC_000071.7",
	"NC_000072.7",
	"NC_000073.7",
	"NC_000074.7",
	"NC_000075.7",
	"NC_000076.7",
	"NC_000077.7",
	"NC_000078.7",
	"NC_000079.7",
	"NC_000080.7",
	"NC_000081.7",
	"NC_000082.7",
	"NC_000083.7",
	"NC_000084.7",
	"NC_000085.7",
	"NC_000086.8",
	"NC_000087.8",
	"NC_005089.1",
	]
```

