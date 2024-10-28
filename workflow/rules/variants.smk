rule call_variants:
    input:
        bam = rules.star_alignreads.output.bam,
        genome = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf
    output:
        vcf = "results/variants/{id}/variants.vcf",
        annotated_vcf = "results/variants/{id}/annotated_variants.vcf"
    params:
        out_dir = "results/variants",
        min_alternate_count = 5,  
        min_coverage = 10,
        snpeff_db = "GRCh38.99" 
    conda:
        "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}.log"
    threads: 8  
    shell: 
        """
        mkdir -p {params.out_dir} && \
        freebayes -f {input.genome} \
            --min-alternate-count {params.min_alternate_count} \
            --min-coverage {params.min_coverage} \
            {input.bam} \
            > {output.vcf} \
            2>> {log} && \
        snpEff ann -v {params.snpeff_db} -Xmx4g {output.vcf} > {output.annotated_vcf} 2>> {log}
        """

rule filter_variants:
    input:
        vcf = rules.call_variants.output.annotated_vcf  
    output:
        vcf_filtered = "results/variants/{id}/20QC_variant.vcf" 
    conda:
        "../envs/bcftools.yml"  
    log:
        "logs/filter_variants_{id}.log"
    shell:
        "bcftools filter -s LowQual -e '%QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"


rule apply_variants:
    input:
        fasta = rules.build_transcriptome.output.transcriptome,  # Remplacez par le chemin réel
        vcf = rules.filter_variants.output.vcf_filtered,  # Fichier VCF filtré
        gtf = rules.download_human_gtf.output.gtf  # Fichier GTF
    output:
        "results/variants/{id}/transcrits_variants.fa"  # Fichier de sortie pour les transcrits avec variantes
    conda:
        "../envs/python.yml"  # Environnement Conda si nécessaire
    log:
        "logs/apply_variants_{id}.log"
    script:
        "../scripts/add_variants.py"  # Chemin vers votre script


# Règle pour l'annotation des variants avec OpenVar
#rule annotate_variants:
#    input:
#        vcf = rules.call_variants.output.vcf
#    output:
#        annotated_vcf = "data/variants/{id}_annotated.vcf"
#    conda:
#        "../envs/openvar.yml"  # Assurez-vous que ce chemin est correct
#    log:
#        "logs/openvar_{id}.log"
#    shell:
#        """
#        openvar -i {input.vcf} -o {output.annotated_vcf} &> {log}
#        """
