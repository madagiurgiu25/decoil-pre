rule decoil:
    input:
        vcf = "{dirname}/sv.sniffles.vcf",
        bw = "{dirname}/coverage.bw",
        bam = bam,
        bai = bam + ".bai",
        ref = ref
    params:
        name = name,
        str = decoilparams,
        dir = "{dirname}"
    output:
        "{dirname}/summary.txt",
        "{dirname}/reconstruct.bed",
        "{dirname}/reconstruct.links.txt",
        "{dirname}/reconstruct.fasta",
        "{dirname}/reconstruct.ecDNA.bed",
        "{dirname}/reconstruct.links.ecDNA.txt",
        "{dirname}/reconstruct.ecDNA.fasta",
        "{dirname}/reconstruct.ecDNA.filtered.bed",
        "{dirname}/reconstruct.links.ecDNA.filtered.txt",
        "{dirname}/reconstruct.ecDNA.filtered.fasta"
    log:
        "{dirname}/logs/logs_decoil"
    shell:
        """
        decoil reconstruct -b {input.bam} -c {input.bw} -i {input.vcf} \
        --name {params.name} -r {input.ref} -o {params.dir} {params.str} &> {log}
        """

        #d.run_reconstruction(input.vcf,
        #                input.bw,
        #                input.bam,
        #                outputdir,
        #                ref,
        #                name=name,
        #                svcaller=svcaller)
