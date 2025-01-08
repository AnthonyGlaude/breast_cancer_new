import pandas as pd
from Bio import SeqIO

def read_vcf(vcf_file):
    mutations = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('##') or line.startswith('#'):
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
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
            mutations.append((chromosome, position, ref_nucleotide, alt_nucleotide, gene_name))
    
    mutations_df = pd.DataFrame(mutations, columns=['chromosome', 'position', 'ref_nucleotide', 'alt_nucleotide', 'gene_name'])
    print("Mutations DataFrame:\n", mutations_df)  # Afficher le DataFrame des mutations
    return mutations_df

def find_exons_for_mutations(vcf_df, gtf_file):
    exons = []
    for _, mutation in vcf_df.iterrows():
        found_exon = False
        with open(gtf_file, 'r') as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] == 'exon':
                    chromosome = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    attributes = fields[8].split(';')
                    gene_name = None

                    for attribute in attributes:
                        if 'gene_name' in attribute:
                            gene_name = attribute.split(' ')[-1].replace('"', '')
                            break

                    if (mutation['chromosome'] == chromosome and 
                        start <= mutation['position'] <= end and 
                        mutation['gene_name'] == gene_name):
                        exons.append({
                            'chromosome': chromosome,
                            'position': mutation['position'],
                            'ref_nucleotide': mutation['ref_nucleotide'],
                            'alt_nucleotide': mutation['alt_nucleotide'],
                            'gene_name': gene_name,
                            'exon_start': start,
                            'exon_end': end
                        })
                        found_exon = True
                        break
        
        if not found_exon:
            print(f"Warning: Mutation at position {mutation['position']} not found in any exon for gene {mutation['gene_name']}.")

    exons_df = pd.DataFrame(exons)
    print("Annotated Exons DataFrame:\n", exons_df)  # Afficher le DataFrame des exons annotés
    return exons_df

def extract_genomic_sequence(fasta_file, exons_df):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    mutated_sequences = []

    for _, exon in exons_df.iterrows():
        chromosome = exon['chromosome']
        position = exon['position']
        alt_nucleotide = exon['alt_nucleotide']
        
        if chromosome in genome:
            seq_record = genome[chromosome]
            exon_sequence = str(seq_record.seq[exon['exon_start'] - 1: exon['exon_end']])
            modified_sequence = exon_sequence[:position - exon['exon_start']] + alt_nucleotide + exon_sequence[position - exon['exon_start'] + 1:]
            
            mutated_sequences.append({
                'chromosome': chromosome,
                'position': position,
                'original_sequence': exon_sequence,
                'modified_sequence': modified_sequence,
                'gene_name': exon['gene_name']
            })
        else:
            print(f"Warning: Chromosome {chromosome} not found in the genome.")

    mutated_sequences_df = pd.DataFrame(mutated_sequences)
    print("Mutated Sequences DataFrame:\n", mutated_sequences_df)  # Afficher le DataFrame des séquences mutées
    return mutated_sequences_df

def modify_transcriptome(transcriptome_file, mutated_exons_df, output_transcriptome_file):
    transcripts = {record.id: str(record.seq) for record in SeqIO.parse(transcriptome_file, 'fasta')}
    
    for _, mutation in mutated_exons_df.iterrows():
        gene_name = mutation['gene_name']
        position = mutation['position']
        alt_nucleotide = mutation['alt_nucleotide']
        
        transcript_id = next((trans_id for trans_id in transcripts if gene_name in trans_id), None)
        
        if transcript_id is not None:
            original_sequence = transcripts[transcript_id]
            transcript_position = position  # À ajuster si nécessaire
            
            if 0 <= transcript_position - 1 < len(original_sequence):
                modified_sequence = (original_sequence[:transcript_position - 1] + 
                                     alt_nucleotide + 
                                     original_sequence[transcript_position:])
                
                new_transcript_id = f"{gene_name}_{mutation['ref_nucleotide']}_to_{mutation['alt_nucleotide']}"
                
                transcripts[transcript_id] = original_sequence  # Séquence normale
                transcripts[new_transcript_id] = modified_sequence  # Séquence modifiée
            else:
                print(f"Warning: Position {transcript_position} is out of range for transcript {transcript_id}.")
        else:
            print(f"Warning: No transcript found for gene {gene_name}.")

    with open(output_transcriptome_file, 'w') as output_handle:
        for trans_id, sequence in transcripts.items():
            output_handle.write(f">{trans_id}\n")
            output_handle.write(f"{sequence}\n")
    print(f"Modified transcriptome saved to {output_transcriptome_file}.")

def validate_and_save_results(mutated_sequences, output_file):
    if mutated_sequences.empty:
        print("No mutated sequences to write to output file.")
        return

    with open(output_file, 'w') as output_fasta:
        for _, seq in mutated_sequences.iterrows():
            gene_name = seq['gene_name']
            modified_sequence = seq['modified_sequence']

            output_fasta.write(f">{gene_name}_modified\n")
            output_fasta.write(f"{modified_sequence}\n")
            
            print(f"Gene: {gene_name}, Modified Sequence: {modified_sequence}")

def main():
    fasta_file = snakemake.input.fasta
    vcf_file = snakemake.input.vcf
    gtf_file = snakemake.input.gtf
    #gtf_file = "/mnt/f/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
    output_transcriptome_file = snakemake.input.transcriptome
    output_file = snakemake.output.fasta

    vcf_df = read_vcf(vcf_file)
    annotated_exons_df = find_exons_for_mutations(vcf_df, gtf_file)
    mutated_sequences_df = extract_genomic_sequence(fasta_file, annotated_exons_df)
    
    modify_transcriptome(fasta_file, annotated_exons_df, output_transcriptome_file)

    validate_and_save_results(mutated_sequences_df, output_file)

if __name__ == "__main__":
    main()
