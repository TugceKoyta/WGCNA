import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os

# === Paramètres ===
base_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/residuals_corrected_lmer/"
covariates_file = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/residuals/sample_covariates.tsv"
conditions_file = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/conditions.txt"
output_dir = os.path.join(base_path, "figures_residuals_clean/")
os.makedirs(output_dir, exist_ok=True)

residual_files = {
            "MG": os.path.join(base_path, "mg_residuals_lmer.parquet"),
            "MT": os.path.join(base_path, "mt_residuals_lmer.parquet"),
                }

palette = {"case": "#E74C3C", "control": "#3498DB"}

# === Lecture du fichier de conditions ===
conditions = pd.read_csv(conditions_file, sep='\t', dtype=str)
conditions.columns = [col.strip().lower() for col in conditions.columns]
conditions = conditions.set_index("sample")

meta = pd.read_csv(covariates_file, sep="\t")
meta = meta.set_index("sample")


# === Traitement pour MG et MT ===
for label, path in residual_files.items():
    print(f"\n== Traitement {label} ==")
    df = pd.read_parquet(path)
    df = df.loc[df.index.intersection(conditions.index)]

    conds = conditions.loc[df.index, "condition"]

    # === Standardisation pour certaines figures ===
    scaled = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)

    # === PCA ===
    pca = PCA(n_components=2)
    components = pca.fit_transform(scaled)
    pca_df = pd.DataFrame(components, columns=["PC1", "PC2"], index=scaled.index)
    pca_df["condition"] = conds

    plt.figure(figsize=(8,6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="condition", palette=palette)
    plt.title(f"PCA - Residuals {label}")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_pca.pdf")
    plt.close()

    # === 🧪 PCA
    family_df = df.copy()  # déjà chargé en haut
    meta_sub = meta.loc[family_df.index]

    pca_family = PCA(n_components=2)
    components_family = pca_family.fit_transform(family_df)

    pca_family_df = pd.DataFrame(components_family, columns=["PC1", "PC2"], index=family_df.index)
    pca_family_df["Family"] = meta_sub["family"].astype(str)

    plt.figure(figsize=(10, 8))
    sns.scatterplot(
            data=pca_family_df,
            x="PC1", y="PC2",
            hue="Family", palette="tab20", s=80
            )

    plt.title(f"PCA - Residuals {label} (Grouped by Family)")
    plt.xlabel(f"PC1 ({pca_family.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca_family.explained_variance_ratio_[1]*100:.1f}%)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Family")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_pca_family.pdf")
    plt.close()
    # === Boxplot top 50 KO ===
    top_ko = df.std().sort_values(ascending=False).head(50).index
    plt.figure(figsize=(10,12))
    sns.boxplot(data=df[top_ko], orient="h", fliersize=1, color="lightgrey")
    plt.title(f"Boxplot - Top 50 KO (residuals) - {label}")
    plt.xlabel("Residuals")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_boxplot_top50.pdf")
    plt.close()

    # === Heatmap ===
    plt.figure(figsize=(14, 10))
    heat_data = scaled[top_ko]
    ax = sns.heatmap(
            heat_data.T,
            cmap="vlag",
            center=0,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": "Z-score"},
            )

    # 🔺 Colorer chaque nom d'échantillon selon la condition
    xticks = ax.get_xticklabels()
    for tick in xticks:
        sample = tick.get_text()
        color = palette.get(conds.loc[sample], "black") if sample in conds.index else "black"
        tick.set_color(color)
        tick.set_rotation(90)

    # 🔻 Options graphiques
    ax.set_title(f"Heatmap - Top 50 KO Residuals - {label}")
    ax.set_xlabel("Samples")
    ax.set_ylabel("KO")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_heatmap_top50.pdf")
    plt.close()


    '''
    row_colors = conds.map(palette)
    plt.figure(figsize=(14, 10))
    sns.heatmap(scaled[top_ko].T, cmap="vlag", center=0,
            xticklabels=True, yticklabels=True, cbar_kws={"label": "Z-score"})
    plt.title(f"Heatmap - Top 50 KO Residuals - {label}")
    plt.xlabel("Samples")
    plt.ylabel("KO")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_heatmap_top50.pdf")
    plt.close()
    '''
    # === Distribution des moyennes ===
    means = df.mean(axis=0)
    plt.figure(figsize=(10,5))
    sns.histplot(means, kde=True, bins=50, color="#2980B9")
    plt.xlim(-50, 200)  # focus sur les KO centraux
    plt.title(f"Distribution of Mean Residuals - {label}")
    plt.xlabel("Average Residuals per KO")
    plt.ylabel("Number of KOs")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{label.lower()}_mean_residuals_hist.pdf")
    plt.close()
