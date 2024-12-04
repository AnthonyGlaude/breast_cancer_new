import pandas as pd
from Bio import SeqIO

def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  
            if line.startswith('#'):
                continue  

            fields = line.strip().split('\t')
            chromosome = fields[0]
            position = int(fields[1])
            ref_nucleotide = fields[3]
            alt_nucleotide = fields[4]
            info_field = fields[7]
            gene_name = '.'  # Valeur par défaut

            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    # Extrait le nom du gène après "ANN=" et le premier '|'
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
            mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name))
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name'])
    return mutations_df

def find_exons_for_mutations(vcf_df, gtf_file):
    exons = []
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':  # On s'intéresse uniquement aux exons
                chromosome = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8].split(';')
                gene_name = None

                # Extraire le nom du gène
                for attribute in attributes:
                    if 'gene_name' in attribute:
                        gene_name = attribute.split(' ')[-1].replace('"', '')
                        break

                # Vérifier si un variant est dans cet exon
                for _, mutation in vcf_df.iterrows():
                    if (mutation['chromosome'] == chromosome and 
                        start <= mutation['position'] <= end and 
                        mutation['gene_name'] == gene_name):  # Vérification du gène
                        exons.append({
                            'chromosome': chromosome,
                            'position': mutation['position'],
                            'ref_nucleotide': mutation['ref_nucleotide'],
                            'alt_nucleotide': mutation['alt_nucleotide'],
                            'gene_name': gene_name,
                            'exon_start': start,
                            'exon_end': end
                        })
                    elif (mutation['chromosome'] == chromosome and 
                        (mutation['position'] < start or mutation['position'] > end) and 
                        mutation['gene_name'] == gene_name):
                        print(f"Warning: Mutation at position {mutation['position']} is outside the exon range ({start}-{end}) for gene {gene_name}.")

    return pd.DataFrame(exons)

def extract_genomic_sequence(fasta_file, exons_df):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    mutated_sequences = []

    for _, exon in exons_df.iterrows():
        chromosome = exon['chromosome']
        position = exon['position']
        alt_nucleotide = exon['alt_nucleotide']
        
        if chromosome in genome:
            seq_record = genome[chromosome]
            # Extraire la séquence de l'exon
            exon_sequence = str(seq_record.seq[exon['exon_start'] - 1: exon['exon_end']])  # -1 pour indexation 0
            # Remplacer le nucléotide à la position de la mutation
            modified_sequence = exon_sequence[:position - exon['exon_start']] + alt_nucleotide + exon_sequence[position - exon['exon_start'] + 1:]
            
            mutated_sequences.append({
                'chromosome': chromosome,
                'position': position,
                'original_sequence': exon_sequence,
                'modified_sequence': modified_sequence,
                'gene_name': exon['gene_name']
            })

    return pd.DataFrame(mutated_sequences)

def modify_transcriptome(transcriptome_file, mutated_exons_df, output_transcriptome_file):
    # Lire le transcriptome
    transcripts = {record.id: str(record.seq) for record in SeqIO.parse(transcriptome_file, 'fasta')}
    
    # Modifier les séquences en fonction des mutations
    for _, mutation in mutated_exons_df.iterrows():
        gene_name = mutation['gene_name']
        position = mutation['position']  # Position dans le génome
        alt_nucleotide = mutation['alt_nucleotide']
        
        # Trouver le transcript correspondant
        transcript_id = next((trans_id for trans_id in transcripts if gene_name in trans_id), None)
        
        if transcript_id is not None:
            original_sequence = transcripts[transcript_id]
            # Convertir la position génomique en position dans le transcript
            # Ceci nécessite une connaissance de la relation entre le génome et le transcript
            transcript_position = position  # Modifier cela si besoin, par exemple, en utilisant un décalage
            
            # Vérifier si la position est valide
            if 0 <= transcript_position - 1 < len(original_sequence):
                # Modifier la séquence
                modified_sequence = (original_sequence[:transcript_position - 1] + 
                                     alt_nucleotide + 
                                     original_sequence[transcript_position:])  # Remplacer le nucléotide
                
                # Renommer le transcript modifié
                new_transcript_id = f"{gene_name}_{mutation['ref_nucleotide']}_to_{mutation['alt_nucleotide']}"
                
                # Ajouter les séquences normales et modifiées au transcriptome
                transcripts[transcript_id] = original_sequence  # Séquence normale
                transcripts[new_transcript_id] = modified_sequence  # Séquence modifiée

    # Écrire le transcriptome modifié dans un fichier FASTA
    with open(output_transcriptome_file, 'w') as output_handle:
        for trans_id, sequence in transcripts.items():
            output_handle.write(f">{trans_id}\n")
            output_handle.write(f"{sequence}\n")

def validate_and_save_results(mutated_sequences, output_file):
    # Écrire les séquences modifiées dans un nouveau fichier FASTA
    with open(output_file, 'w') as output_fasta:
        for seq in mutated_sequences:
            gene_name = seq['gene_name']
            modified_sequence = seq['modified_sequence']

            # Enregistrer la séquence modifiée avec des annotations
            output_fasta.write(f">{gene_name}_modified\n")
            output_fasta.write(f"{modified_sequence}\n")
            
            # Vérification de chaque modification
            print(f"Gene: {gene_name}, Modified Sequence: {modified_sequence}")

def main():
    fasta_file = snakemake.input.fasta
    vcf_file = snakemake.input.vcf
    gtf_file = snakemake.input.gtf
    output_file = snakemake.output.fasta
    output_transcriptome_file = snakemake.output.transcriptome

    vcf_df = read_vcf(vcf_file)
    annotated_exons_df = find_exons_for_mutations(vcf_df, gtf_file)
    mutated_sequences_df = extract_genomic_sequence(fasta_file, annotated_exons_df)
    
    modify_transcriptome(fasta_file, annotated_exons_df, output_transcriptome_file)

    # Validation et enregistrement des résultats
    validate_and_save_results(mutated_sequences_df, output_file)

if __name__ == "__main__":
    main()
