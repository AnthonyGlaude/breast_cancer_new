import os
import pandas as pd
import matplotlib.pyplot as plt

# Chemins des fichiers
vcf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, low_memory=False)

gtf_df.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NEW_COLUMN']

def get_feature_for_variant(row):
    chrom = row['#CHROM']
    pos = row['POS']

    # Filtrer le DataFrame GTF pour trouver les fonctionnalités correspondantes
    matching_gtf = gtf_df[(gtf_df['CHROM'] == chrom) & (gtf_df['start'] <= pos) & (gtf_df['end'] >= pos)]

    if not matching_gtf.empty:
        feature_counts = matching_gtf['feature'].value_counts()
        
        # Liste des fonctionnalités prioritaires
        priority_features = ['transcript', '5-utr', '3-utr', 'exon', 'gene']

        # Vérifier les fonctionnalités en fonction de l'ordre de priorité
        for feature in priority_features:
            if feature in feature_counts.index:
                return feature  # Retourne la première fonctionnalité prioritaire rencontrée

        # Gérer les fonctionnalités non prioritaires
        all_features = feature_counts.index.tolist()
        other_features = [feature for feature in all_features if feature not in priority_features]

        if other_features:
            return other_features  # Retourne les fonctionnalités autres trouvées

        # Si aucune fonctionnalité prioritaire ou autre n'est trouvée, retourner la plus fréquente
        return feature_counts.idxmax() if not feature_counts.empty else "Aucune fonctionnalité"

    return "Région Intragénique"  # Si aucune correspondance

def main():
    global gtf_df
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, low_memory=False)
    gtf_df.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    variants_folder = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/results/variants/"
    total_counts = {"transcript": 0, "exon": 0, "gene": 0, "5-utr": 0, "3-utr": 0, "Région Intragénique": 0}

    for patient_folder in os.listdir(variants_folder):
        patient_path = os.path.join(variants_folder, patient_folder)
        vcf_file = os.path.join(patient_path, "20QC_variant.vcf")

        if os.path.isfile(vcf_file):
            try:
                vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, low_memory=False, on_bad_lines='skip')

                for index, row in vcf_df.iterrows():
                    print(f"Traitement du variant: {row}")  # Impression pour le débogage
                    feature = get_feature_for_variant(row)
                    print(f"Fonctionnalité trouvée: {feature}")  # Impression pour le débogage

                    if isinstance(feature, list):
                        for f in feature:
                            total_counts[f] = total_counts.get(f, 0) + 1
                    else:
                        total_counts[feature] += 1
            except Exception as e:
                print(f"Erreur lors de la lecture du fichier {vcf_file}: {e}")

    print("Comptes totaux:", total_counts)  # Debugging

    plt.figure(figsize=(8, 6), dpi=600)
    plt.pie(total_counts.values(), labels=total_counts.keys(), autopct='%1.1f%%', startangle=90)
    plt.ylabel('')
    plt.title(f"Distribution des variants par localisation intragénique (N={sum(total_counts.values())})")
    plt.savefig("piechart_high_res.png", dpi=600)
    plt.show()

if __name__ == "__main__":
    main()
