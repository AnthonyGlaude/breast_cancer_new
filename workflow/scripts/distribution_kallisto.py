import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np  # Importation de numpy nécessaire pour les calculs trigonométriques

# Fonction pour charger les données VCF avec pandas
def load_vcf_with_pandas(vcf_file):
    """
    Charge les données VCF à partir du fichier en utilisant pandas
    """
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [f'Sample{i+1}' for i in range(vcf_df.shape[1] - 9)]
    return vcf_df

# Fonction pour extraire le type de variant de la colonne 'INFO'
def extract_variant_type(info_column):
    """
    Extrait le type de variant à partir de la colonne 'INFO' du fichier VCF.
    """
    info_str = info_column
    type_info = [entry.split('=')[1] for entry in info_str.split(';') if entry.startswith('TYPE=')]
    return type_info[0] if type_info else 'Autre'

# Fonction pour afficher les types de variants détectés avec une légende
def plot_variant_types_with_legend(variants_df, output_plot_path):
    """
    Crée un pie chart des différents types de variants détectés avec une légende colorée.
    """
    # Extraire les types de variants de la colonne 'INFO'
    variants_df['variant_type'] = variants_df['INFO'].apply(extract_variant_type)
    
    # Compter les types de variants
    variant_types = variants_df['variant_type'].value_counts()
    
    # Créer un pie chart sans pourcentage ou labels
    fig, ax = plt.subplots(figsize=(7, 7))  # Créer une figure plus grande pour éviter l'encombrement
    wedges, _ = ax.pie(variant_types, 
                       startangle=90, 
                       colors=sns.color_palette("Set1", len(variant_types)),
                       wedgeprops={'edgecolor': 'black', 'linewidth': 0.5})
    
    # Créer une légende avec les couleurs et les labels (types de variants et leurs pourcentages)
    legend_labels = [f"{variant_types.index[i]}: {variant_types.iloc[i]} ({variant_types.iloc[i] / variant_types.sum() * 100:.1f}%)" 
                     for i in range(len(variant_types))]

    ax.legend(wedges, legend_labels, title="Types de Variants", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=12, frameon=False)  # Ajouter la légende à droite du graphique

    # Ajouter un titre
    ax.set_title("Types de Variants Détectés", fontsize=16)
    
    # Sauvegarder le graphique avec un DPI élevé
    plt.savefig(output_plot_path, dpi=720)  # Augmenter la résolution à 300 DPI
    plt.close()

# Fonction pour générer un graphique de distribution des abondances
def plot_abundance_distribution(abundance_data, output_plot_path):
    """
    Crée un graphique de distribution des abondances.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(abundance_data, kde=True, color='blue', bins=30)
    plt.title("Distribution des Abondances des Variants", fontsize=16)
    plt.xlabel("Abondance", fontsize=12)
    plt.ylabel("Fréquence", fontsize=12)
    
    # Sauvegarder le graphique avec un DPI élevé
    plt.savefig(output_plot_path, dpi=720)  # Augmenter la résolution à 300 DPI
    plt.close()

# Fonction pour afficher les types de variants avec un pie chart détaillé
def plot_variant_types_piechart(variants_df, output_plot_path):
    """
    Crée un pie chart des types de variants détectés, incluant les labels avec ajustement.
    """
    # Extraire les types de variants de la colonne 'INFO'
    variants_df['variant_type'] = variants_df['INFO'].apply(extract_variant_type)
    
    # Compter les types de variants
    variant_types = variants_df['variant_type'].value_counts()
    
    # Créer un pie chart avec des labels ajustés
    fig, ax = plt.subplots(figsize=(8, 8))  # Agrandir la taille de la figure
    wedges, texts, autotexts = ax.pie(variant_types, 
                                      labels=variant_types.index,
                                      autopct='%1.1f%%',
                                      startangle=90, 
                                      colors=sns.color_palette("Set1", len(variant_types)),
                                      wedgeprops={'edgecolor': 'black', 'linewidth': 0.5})
    
    # Personnaliser les tailles des labels et textes pour éviter le chevauchement
    for text in texts:
        text.set_fontsize(10)
        text.set_horizontalalignment('center')
    
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_horizontalalignment('center')
    
    # Ajouter un titre
    ax.set_title("Types de Variants Détectés", fontsize=16)
    
    # Ajouter la légende à droite du graphique
    ax.legend(wedges, [f"{variant_types.index[i]}: {variant_types.iloc[i]} ({variant_types.iloc[i] / variant_types.sum() * 100:.1f}%)" 
                       for i in range(len(variant_types))],
              title="Types de Variants", loc="center left", bbox_to_anchor=(1.05, 0.5),
              fontsize=12, frameon=False)  # Placer la légende à droite du graphique

    # Optimiser l'espace autour du graphique avec tight_layout()
    plt.tight_layout()

    # Sauvegarder le graphique avec un DPI élevé
    plt.savefig(output_plot_path, dpi=720)  # Augmenter la résolution à 300 DPI
    plt.close()

# Exemple pour les variants (en supposant un fichier de variants '20QC_variant.vcf')
variants_file = "F:/breast_cancer/workflow/results/fraction/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"

# Charger les données VCF avec pandas
vcf_df = load_vcf_with_pandas(variants_file)

# Générer le graphique pour les types de variants détectés avec la légende
output_plot_path_variant_types_legend = "F:/breast_cancer/workflow/results/fraction/dge/kallisto/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/variant_types_with_legend.png"
plot_variant_types_with_legend(vcf_df, output_plot_path_variant_types_legend)

# Générer le graphique pour la distribution des abondances
abundance_data = np.random.uniform(0, 100, size=1000)  # Exemple de données d'abondance
output_plot_path_abundance = "F:/breast_cancer/workflow/results/fraction/dge/kallisto/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/abundance_distribution.png"
plot_abundance_distribution(abundance_data, output_plot_path_abundance)

# Générer le graphique pour les types de variants en piechart détaillé
output_plot_path_variant_types_piechart = "F:/breast_cancer/workflow/results/fraction/dge/kallisto/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/variant_types_piechart.png"
plot_variant_types_piechart(vcf_df, output_plot_path_variant_types_piechart)

print(f"Graphiques générés et enregistrés sous :")
print(f"- Types de variants avec légende : {output_plot_path_variant_types_legend}")
print(f"- Distribution des abondances : {output_plot_path_abundance}")
print(f"- Types de variants (pie chart détaillé) : {output_plot_path_variant_types_piechart}")
