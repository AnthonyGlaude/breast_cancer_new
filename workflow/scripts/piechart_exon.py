import pandas as pd
import matplotlib.pyplot as plt

# Chemins des fichiers
gtf_file = "F:/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
vcf_file = "F:/breast_cancer/workflow/results/fraction/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)

gtf_df.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NEW_COLUMN']

# Fonction pour obtenir les attributs comme 'gene_id' et 'transcript_id'
def extract_attributes(attribute_str):
    attributes = {}
    for attr in attribute_str.split(';'):
        if attr.strip():  # Ignore les éléments vides
            key, value = attr.split(' ', 1)
            key = key.strip().lower()  # Assure-toi que la clé est en minuscules
            value = value.strip().strip('"')  # Retirer les guillemets autour de la valeur
            attributes[key] = value
    return attributes

# Applique l'extraction des attributs sur la colonne 'attribute'
gtf_df['attributes'] = gtf_df['attribute'].apply(extract_attributes)

# Extraire 'gene_id' et 'transcript_id'
gtf_df['gene_id'] = gtf_df['attributes'].apply(lambda x: x.get('gene_id', None))
gtf_df['transcript_id'] = gtf_df['attributes'].apply(lambda x: x.get('transcript_id', None))

# Fonction pour obtenir la région intragénique pour chaque variant
def get_feature_for_variant(row):
    chrom = row['#CHROM']
    pos = row['POS']
    
    # Debugging: Afficher les valeurs de chrom et pos pour chaque variant
    print(f"Processing variant at Chrom: {chrom}, Position: {pos}")
    
    # Rechercher les régions dans le GTF qui chevauchent la position du variant
    matching_gtf = gtf_df[(gtf_df['CHROM'] == chrom) & (gtf_df['start'] <= pos) & (gtf_df['end'] >= pos)]
    
    # Debugging: Afficher le résultat de la recherche
    if not matching_gtf.empty:
        print(f"Match found: {matching_gtf[['feature', 'start', 'end']]}")

        # Comptage des occurrences de chaque fonctionnalité
        feature_counts = matching_gtf['feature'].value_counts()

        # Si transcript, exon et gene sont présents, appliquer la règle : transcript > exon > gene
        if 'transcript' in feature_counts.index:
            return "transcript"
        elif 'exon' in feature_counts.index:
            return "exon"
        elif 'gene' in feature_counts.index:
            return "gene"
        else:
            return feature_counts.idxmax()  # Si aucun des précédents, retourner la fonctionnalité la plus fréquente
    else:
        # Debugging: Afficher quand aucune correspondance n'est trouvée
        print(f"No match found for variant at Chrom: {chrom}, Position: {pos}")
    
    return "Région Intragénique"  # Si aucune correspondance

# Applique la fonction sur chaque ligne du DataFrame VCF
vcf_df['Feature'] = vcf_df.apply(get_feature_for_variant, axis=1)
feature_counts = vcf_df['Feature'].value_counts()

# Définir les couleurs spécifiques pour chaque catégorie
colors = {
    "transcript": "#548235",      # Couleur pour transcript
    "exon": "#FFD700",            # Couleur pour exon (exemple)
    "gene": "#98FB98",            # Couleur pour gene (exemple)
    "Région Intragénique": "#A3D08E"  # Couleur pour Région Intragénique
}

# Plot pie chart avec couleurs spécifiques
plt.figure(figsize=(8, 6), dpi=600)
feature_counts.plot(kind='pie', autopct='%1.1f%%', startangle=90, cmap='Set3', colors=[colors.get(label, "#D3D3D3") for label in feature_counts.index], legend=False)

plt.ylabel('')
#plt.title("Distribution des variants par localisation intragénique")
plt.show()
plt.savefig("piechart_high_res.png", dpi=600)

# Affiche les comptages des différentes régions
print(feature_counts)
