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
