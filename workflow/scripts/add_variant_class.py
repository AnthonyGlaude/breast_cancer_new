from Bio import SeqIO
import pandas as pd


class Mutation:
    def __init__(self, chromosome, position, ref_nucleotide, alt_nucleotide, gene_name='.'):
        self.chromosome = chromosome
        self.position = position
        self.ref_nucleotide = ref_nucleotide
        self.alt_nucleotide = alt_nucleotide
        self.gene_name = gene_name

    def __repr__(self):
        return (f"Mutation(chromosome={self.chromosome}, position={self.position}, "
                f"ref={self.ref_nucleotide}, alt={self.alt_nucleotide}, gene={self.gene_name})")


class Exon:
    def __init__(self, chromosome, start, end, gene_name):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.gene_name = gene_name

    def contains(self, mutation):
        """Check if a mutation is within this exon."""
        return (self.chromosome == mutation.chromosome and
                self.start <= mutation.position <= self.end and
                self.gene_name == mutation.gene_name)

    def __repr__(self):
        return (f"Exon(chromosome={self.chromosome}, start={self.start}, "
                f"end={self.end}, gene={self.gene_name})")


class Genome:
    def __init__(self, fasta_file):
        self.genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    def get_sequence(self, chromosome, start, end):
        """Extract a sequence from the genome."""
        if chromosome in self.genome:
            seq_record = self.genome[chromosome]
            return str(seq_record.seq[start - 1:end])
        else:
            raise ValueError(f"Chromosome {chromosome} not found in the genome.")


class Annotator:
    def __init__(self, gtf_file):
        self.exons = self._parse_gtf(gtf_file)

    def _parse_gtf(self, gtf_file):
        exons = []
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

                    exons.append(Exon(chromosome, start, end, gene_name))
        return exons

    def find_exons_for_mutation(self, mutation):
        """Find the exon that contains the given mutation."""
        for exon in self.exons:
            if exon.contains(mutation):
                return exon
        return None


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
            gene_name = '.'  # Default value

            for item in info_field.split(';'):
                if item.startswith('ANN='):
                    gene_info = item.split('|')
                    if len(gene_info) > 3:
                        gene_name = gene_info[3]
            mutations.append(Mutation(chromosome, position, ref_nucleotide, alt_nucleotide, gene_name))

    return mutations


def annotate_mutations(vcf_file, gtf_file):
    mutations = read_vcf(vcf_file)
    annotator = Annotator(gtf_file)

    annotated_mutations = []
    for mutation in mutations:
        exon = annotator.find_exons_for_mutation(mutation)
        if exon:
            annotated_mutations.append({
                'chromosome': mutation.chromosome,
                'position': mutation.position,
                'ref_nucleotide': mutation.ref_nucleotide,
                'alt_nucleotide': mutation.alt_nucleotide,
                'gene_name': mutation.gene_name,
                'exon_start': exon.start,
                'exon_end': exon.end
            })
        else:
            print(f"Warning: Mutation {mutation} not found in any exon.")

    return pd.DataFrame(annotated_mutations)


def extract_modified_sequences(fasta_file, annotated_mutations_df):
    genome = Genome(fasta_file)
    mutated_sequences = []

    for _, row in annotated_mutations_df.iterrows():
        chromosome = row['chromosome']
        position = row['position']
        ref_nucleotide = row['ref_nucleotide']
        alt_nucleotide = row['alt_nucleotide']
        exon_start = row['exon_start']
        exon_end = row['exon_end']
        gene_name = row['gene_name']

        try:
            exon_sequence = genome.get_sequence(chromosome, exon_start, exon_end)
            modified_sequence = (exon_sequence[:position - exon_start] + 
                                 alt_nucleotide + 
                                 exon_sequence[position - exon_start + 1:])
            mutated_sequences.append({
                'chromosome': chromosome,
                'position': position,
                'original_sequence': exon_sequence,
                'modified_sequence': modified_sequence,
                'gene_name': gene_name
            })
        except ValueError as e:
            print(e)

    return pd.DataFrame(mutated_sequences)


def main():
    fasta_file = snakemake.input.fasta
    vcf_file = snakemake.input.vcf
    gtf_file = snakemake.input.gtf
    output_file = snakemake.output.fasta

    annotated_mutations_df = annotate_mutations(vcf_file, gtf_file)
    mutated_sequences_df = extract_modified_sequences(fasta_file, annotated_mutations_df)
    mutated_sequences_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
