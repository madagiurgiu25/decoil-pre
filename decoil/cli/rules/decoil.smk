rule decoil:
    input:
        vcf = "{dirname}/sv.sniffles.vcf",
        bw = "{dirname}/coverage.bw",
        bam = bam,
        bai = bam + ".bai",
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
    run:
        d.run_reconstruction(input.vcf,
                        input.bw,
                        input.bam,
                        outputdir,
                        ref,
                        name=name,
                        svcaller=svcaller)