rule svcalling:
        input:
            bam
        output:
            "{dirname}/sv.sniffles.vcf"
        log:
            "{dirname}/logs/logs_sniffles"
        params:
            threads = threads
        shell:
            """
            sniffles -t {params.threads} -m {input} \
            -v {output} --min_homo_af 0.7 \
            --min_het_af 0.1 --min_length 50 \
            --cluster --genotype --min_support 4 --report-seq &> {log}
            """

            # """
            # sniffles -t {params.threads} -i {bam} -v {output} --minsvlen 50 --minsupport 4 &> {log}
            # """

rule survivor:
    input:
        "{dirname}/sv.sniffles.vcf"
    output:
        "{dirname}/sv.sniffles.bedpe"
    log:
        "{dirname}/logs/logs_survivor"
    shell:
        """
        SURVIVOR vcftobed {input} -1 -1 {output} &> {log}
        """

rule coverage:
    input:
        bam=bam,
        bai=bam + ".bai"
    output:
        "{dirname}/coverage.bw"
    # conda:
    #     "envs/cov.yaml"
    log:
        "{dirname}/logs/logs_cov"
    params:
        threads=threads
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 -p {params.threads} &> {log}
        """
