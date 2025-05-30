rule download_human_gtf:
    """ Download gtf of human genome from Ensembl """
    output:
        gtf = 'data/references/gtf/Homo_sapiens.GRCh38.110.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        """
        mkdir -p data/references/gtf &&
        wget -O temp_gtf.gz {params.link} &&
        gunzip temp_gtf.gz &&
        mv temp_gtf {output.gtf} &&  # Renommer le fichier temporaire en gtf
        rm -f temp_gtf.gz
        """


rule download_human_genome:
    """Download the reference genome (fasta file) of human from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/homo_sapiens_genome.fa'
    params:
        link = config['download']['human_genome_fa']
    shell:
        """
        mkdir -p data/references/genome_fa && 
        wget -O temp_genome.gz {params.link} &&
        gunzip temp_genome.gz && 
        mv temp_genome {output.genome} &&
        rm -f temp_genome.gz
        """


#rule download_human_gff3:
#   """ Download gff3 of human genome from Ensembl """
#    output:
#        gff3 = 'data/references/gff3/homo_sapiens.gff3'
#    params:
#        link = config['download']['human_gff3']
#    shell:
#        """
#        mkdir -p data/references/gff3 && 
#        wget -O temp.gz {params.link} && 
#        gunzip temp.gz && 
#        mv temp {output.gff3} &&
#        rm temp.gz
#        """

