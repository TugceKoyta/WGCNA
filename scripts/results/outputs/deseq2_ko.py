import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from bioservices import KEGG
import os

# === 📁 CHEMINS D'ACCÈS ===
base_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
counts_files = {
            "mg": os.path.join(base_path, "mg_ko_raw_counts_matrix.parquet"),
            "mt": os.path.join(base_path, "mt_ko_raw_counts_matrix.parquet")
                }
conditions_file = os.path.join(base_path, "conditions.txt")
output_dir = os.path.join(base_path, "deseq2_results2/")
os.makedirs(output_dir, exist_ok=True)

# === ⚙️ PARAMÈTRES ===
log2fc_thresh = 0.5
#padj_thresh = 0.1
padj_thresh_dict = {
        "mg": 0.9,
        "mt": 0.1  # par exemple
        }

# === 🔁 TRAITEMENT POUR MG / MT ===
for label, path in counts_files.items():
    print(f"\n🔬 Traitement DESeq2 pour {label.upper()}")
    df = pd.read_parquet(path)
    print("{label} matrix:", df.shape)
    print("{label} - KO with non-zero sum:", (df.sum(axis=1) > 0).sum())

    if df.index.name != "Geneid":
        df.index.name = "Geneid"

    df_ko = df
    df_ko.index = df_ko.index.astype(str)

    #df_ko = df_ko.T  # KO en lignes, échantillons en colonnes

    # 📂 Chargement des conditions
    conds = pd.read_csv(conditions_file, sep=r'\s+', header=0)
    conds["sample"] = conds["sample"].astype(str).str.strip()
    conds = conds.set_index("sample")

    df_ko = df_ko.loc[:, df_ko.columns.isin(conds.index)]
    conds = conds.loc[df_ko.columns]
    print(f"✅ {df_ko.shape[1]} échantillons, {df_ko.shape[0]} KO")


    # Nettoyage des noms d’échantillons
    df_ko.columns = df_ko.columns.astype(str).str.strip().str.replace('\r', '').str.replace('\t', '')
    conds.index = conds.index.astype(str).str.strip().str.replace('\r', '').str.replace('\t', '')

    # Diagnostic
    print("Échantillons dans la matrice :", df_ko.columns.tolist()[:5])
    print("Échantillons dans conditions.txt :", conds.index.tolist()[:5])

    # Intersection
    common_samples = df_ko.columns.intersection(conds.index)
    print(f"🎯 {len(common_samples)} échantillons en commun.")

    # Filtrage
    if len(common_samples) < 2:
        print(f"❌ Moins de 2 échantillons avec conditions valides pour {label.upper()} -> ignoré")
        continue

    df_ko = df_ko[common_samples]
    conds = conds.loc[common_samples]


    # 🔍 Filtrer les échantillons communs
    common_samples = df_ko.columns.intersection(conds.index)
    if len(common_samples) < 2:
        print(f"❌ Moins de 2 échantillons avec conditions valides pour {label.upper()} -> ignoré")
        continue

    df_ko = df_ko[common_samples]
    conds = conds.loc[common_samples]

    if conds['condition'].nunique() < 2:
        print(f"❌ Moins de deux conditions uniques dans les métadonnées pour {label.upper()} -> ignoré")
        continue

    print(f"✅ {df_ko.shape[1]} échantillons, {df_ko.shape[0]} KO")

    # 🧬 DESeq2
    dds = DeseqDataSet(
            counts=df_ko.T.astype(int),
            metadata=conds,
            design_factors=["condition"]
            )

    dds.deseq2()
    padj_thresh = padj_thresh_dict[label]

    stats = DeseqStats(dds, alpha=padj_thresh)
    stats.summary()
    res_df = stats.results_df.copy().dropna(subset=["padj"])
    res_df_filtered = res_df[(res_df["padj"] < padj_thresh) & (abs(res_df["log2FoldChange"]) > log2fc_thresh)]
    res_df.to_csv(os.path.join(output_dir, f"{label}_ko_deseq2_results.csv"))
    res_df_filtered.to_csv(os.path.join(output_dir, f"{label}_ko_significatifs.csv"))
    print(f"✅ {res_df_filtered.shape[0]} KO significatifs ({label.upper()})")

    # 📈 Volcano Plot
    res_df["-log10(padj)"] = -np.log10(res_df["padj"])
    plt.figure(figsize=(10, 6))
    plt.scatter(res_df['log2FoldChange'], res_df['-log10(padj)'], alpha=0.3, color='grey')
    plt.scatter(res_df_filtered['log2FoldChange'], -np.log10(res_df_filtered['padj']), color='red', alpha=0.6)
    plt.axhline(-np.log10(padj_thresh), color='blue', linestyle='--')
    plt.axvline(-log2fc_thresh, color='green', linestyle='--')
    plt.axvline(log2fc_thresh, color='green', linestyle='--')
    plt.title(f"Volcano plot - {label.upper()}")
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10 Adjusted p-value")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{label}_volcano_plot.png"), dpi=300)
    plt.close()

    # MA Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(res_df["baseMean"], res_df["log2FoldChange"], alpha=0.3, color="grey")
    plt.scatter(res_df_filtered["baseMean"], res_df_filtered["log2FoldChange"], color="red")
    plt.xscale("log")
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel("baseMean (log scale)")
    plt.ylabel("log2 Fold Change")
    plt.title(f"MA Plot - {label.upper()}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ma_plot.png"), dpi=300)
    plt.close()
