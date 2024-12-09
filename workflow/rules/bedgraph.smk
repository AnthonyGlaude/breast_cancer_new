rule vcf_to_bedgraph:
    input:
        vcf="/mnt/f/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf",
        genome=rules.download_human_genome.output.genome
    output:
        bedgraph="results/graph/{id}.bedgraph"
    conda:
        "../envs/vcf_to_bedgraph.yml"
    shell:
        """
        # Index du génome FASTA (si non déjà indexé)
        if [ ! -f {input.genome}.fai ]; then
            samtools faidx {input.genome}
        fi

        # Conversion VCF -> BED (positions des variants)
        bcftools query -f '%CHROM\t%POS\t%END\t%QUAL\n' {input.vcf} > results/{wildcards.id}_variants.bed

        # Création du fichier BEDGraph
        bedtools genomecov -bg -i results/{wildcards.id}_variants.bed -g {input.genome}.fai > {output.bedgraph}
        """

rule create_chrom_sizes:
    input:
        gtf = rules.download_human_gtf.output.gtf
    output:
        chrom_sizes = "results/graph/chrom.sizes"
    shell:
        """
        awk '$3 == "gene" {{print $1, $5}}' {input.gtf} |
        sort -k1,1 -k2,2n | 
        uniq |
        awk '{{print $1, $2}}' > {output.chrom_sizes}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "results/graph/{id}.bedgraph",
        chrom_sizes = rules.create_chrom_sizes.output.chrom_sizes
    output:
        bw = "results/graph/{id}.bw"
    conda:
        "../envs/vcf_to_bedgraph.yml"
    shell:
        """
        bedGraphToBigWig {input.bedgraph} {input.chrom_sizes} {output.bw}
        """
