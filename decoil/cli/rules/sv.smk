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

if svcaller == "sniffles1":

    rule svcalling:
        input:
            bam
        output:
            "{dirname}/sv.sniffles.vcf"
        #conda:
        #    get_conda_env_sv()
        log:
            "{dirname}/logs/logs_sniffles1"
        params:
            threads = threads
        shell:
            """
            sniffles -t {params.threads} -m {input} \
            -v {output} --min_homo_af 0.7 \
            --min_het_af 0.1 --min_length 50 \
            --cluster --genotype --min_support 4 --report-seq &> {log}
            """

elif svcaller == "sniffles2":

    rule svcalling:
        input:
            bam
        output:
            "{dirname}/sv.sniffles.vcf"
        conda:
            "../envs/sniffles2.yaml"
        log:
            "{dirname}/logs/logs_sniffles2"
        params:
            threads = threads
        shell:
            """
            sniffles --threads {params.threads} --input {input} \
            --vcf {output} --no-qc --minsupport 4 --no-qc --mosaic \
            --mosaic-include-germline --mosaic-af-max 0.4 --mosaic-af-min 0.01 --allow-overwrite &> {log}
            """

rule survivor:
    input:
        "{dirname}/sv.sniffles.vcf"
    output:
        "{dirname}/sv.sniffles.bedpe"
    conda:
            "../envs/survivor.yaml"
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
