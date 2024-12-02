import platform

# Function to select the environment file based on the system
def get_conda_env_sv():
    system = platform.system()
    if system == "Linux":
        return "../envs/sv.yaml"
    elif system == "Darwin":  # macOS
        return "../envs/sv_macos.yaml"
    elif system == "Windows":
        return "../envs/sv.yaml"
    else:
        raise ValueError("Unsupported system: " + system)

rule svcalling:
    input:
        bam = bam,
        ref = ref
    output:
        "{dirname}/sv.sniffles.vcf"
#    conda:
#        get_conda_env_sv()
    log:
        "{dirname}/logs/logs_sniffles"
    params:
        threads = threads,
        svcaller = svcaller
    run:
        if svcaller == "sniffles1":
            shell("""sniffles -t {params.threads} -m {input.bam} \
            -v {output} --min_homo_af 0.7 \
            --min_het_af 0.1 --min_length 50 \
            --cluster --genotype --min_support 4 --report-seq &> {log}""")

        elif svcaller == "sniffles2" and filt_version == 1:
            shell("""sniffles --threads {params.threads} --input {input.bam} --qc-coverage 4 \
            --reference {input.ref} --vcf {output} --no-qc --minsupport 4 --allow-overwrite &> {log}""")

        elif svcaller == "sniffles2" and filt_version == 2:
            shell("""sniffles --threads {params.threads} --input {input.bam} --qc-coverage 4 \
            --reference {input.ref} --vcf {output} --no-qc --minsupport 4 --mosaic --minsvlen 500 \
            --mosaic-include-germline --mosaic-af-max 0.4 --mosaic-af-min 0.01 --allow-overwrite &> {log}""")

        elif svcaller == "sniffles2" and filt_version == 3:
            shell("""sniffles --threads {params.threads} --input {input.bam} --qc-coverage 4 \
            --reference {input.ref} --minsvlen 500 --cluster-merge-pos 300 --vcf {output} \
            --allow-overwrite --phase &> {log}
            """)
        else:
            shell("""> {output}""")
      

rule survivor:
    input:
        "{dirname}/sv.sniffles.vcf"
    output:
        "{dirname}/sv.sniffles.bedpe"
#    conda:
#        "../envs/survivor.yaml"
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
    log:
        "{dirname}/logs/logs_cov"
    params:
        threads=threads
    shell:
        """
        bamCoverage --bam {input.bam} -o {output} --binSize 50 -p {params.threads} &> {log}
        """
