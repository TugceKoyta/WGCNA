import pandas as pd
import numpy as np
import os
import gzip
from collections import defaultdict
import glob

# === 📁 CHEMINS D'ACCÈS ===
base_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
counts_files = {
            "mg": os.path.join(base_path, "mg_feature_table.parquet"),
            "mt": os.path.join(base_path, "mt_feature_table.parquet")
                }
gff_files = glob.glob("/scratch/users/tkoytaviloglu/data/raw/MUST/M*/Analysis/annotation/annotation_CDS_RNA_hmms.gff.gz"
)
output_files = {
            "mg": os.path.join(base_path, "ko_mg_gff_tpm.parquet"),
            "mt": os.path.join(base_path, "ko_mt_gff_tpm.parquet")
                }

# === ⚙️ PARAMÈTRES ===
#abundance_threshold = 30
#presence_threshold = 3
presence_fraction = 0.2  # 80% des échantillons

# === 🔧 FONCTIONS ===
def parse_gff(gff_files):
    #Extrait la longueur totale des gènes associés à des KO (kegg_ko=Kxxxxx) à partir d'une liste de fichiers GFF.
    gene_ko_lengths = defaultdict(list)
    for gff_file in gff_files:
        with gzip.open(gff_file, 'rt') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                start = int(parts[3])
                end = int(parts[4])
                attributes = parts[8]
                ko_id = None
                for attr in attributes.split(';'):
                    if attr.strip().startswith("kegg_ko="):
                        ko_id = attr.strip().split('=')[1]
                        break
                if ko_id:
                    gene_ko_lengths[ko_id].append((start, end))

    # Calcul des longueurs non chevauchantes pour chaque KO
    ko_lengths = {}
    for ko_id, exons in gene_ko_lengths.items():
        sorted_exons = sorted(exons, key=lambda x: x[0])
        merged = []
        current_start, current_end = sorted_exons[0]
        for start, end in sorted_exons[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
                merged.append((current_start, current_end))
                total_length = sum(end - start + 1 for start, end in merged)
                ko_lengths[ko_id] = total_length
    return ko_lengths


def calculate_tpm(counts_df, gene_lengths):
    valid = counts_df.index.intersection(gene_lengths.keys())
    counts_df = counts_df.loc[valid]
    lengths = pd.Series(gene_lengths).loc[valid]
    rpk = counts_df.div(lengths, axis=0) * 1e3
    scaling = rpk.sum(axis=0) / 1e6
    tpm = rpk.div(scaling, axis=1)
    return tpm

#def filter_tpm(tpm_df, abundance_threshold, presence_threshold, presence_fraction):
    #abundance_filter = (tpm_df >= abundance_threshold).sum(axis=1) >= presence_threshold
    #filtered = tpm_df.loc[abundance_filter]
    #presence_filter = (filtered > 0).sum(axis=1) >= presence_fraction * filtered.shape[1]
    #return filtered.loc[presence_filter]
def filter_tpm(tpm_df, presence_fraction):
    presence_filter = (tpm_df > 0).sum(axis=1) >= presence_fraction * tpm_df.shape[1]
    return tpm_df.loc[presence_filter]


# === 🔁 TRAITEMENT POUR MG / MT ===
for label in ["mg", "mt"]:
    print(f"\n🔬 Traitement des données {label.upper()}")



    df = pd.read_parquet(counts_files[label])
    if df.index.name != "index":
        if "index" in df.columns:
            df = df.set_index("index")
    ko_cols = [col for col in df.columns if col.startswith(f"{label}_mantis_kegg_ko_counts_")]
    if not ko_cols:
        print(f"❌ Aucune colonne KO détectée pour {label.upper()}")
        continue

    df_ko = df[ko_cols]
    df_ko.columns = df_ko.columns.str.replace(f"{label}_mantis_kegg_ko_counts_", "", regex=False)
    df_ko.index = df_ko.index.astype(str)

    print("📥 Extraction des longueurs de gènes...")
    gene_lengths = parse_gff(gff_files)
    print("🧮 Calcul des TPM...")
    tpm_df = calculate_tpm(df_ko.T, gene_lengths)
    print("🔍 Application des filtres...")
    # === 🧼 Nettoyage complémentaire des KO avant filtrage final ===

    # 1. Filtrage par médiane (KO peu exprimés dans la majorité des échantillons)
    median_filter = tpm_df.median(axis=1) > 0
    tpm_df = tpm_df.loc[median_filter]

    # 2. Suppression des KO avec des valeurs extrêmement élevées dans quelques échantillons
    extreme_threshold = 10000  # seuil TPM considéré comme extrême
    max_extreme_occurence = 10

    extreme_occurrences = (tpm_df > extreme_threshold).sum(axis=1)
    tpm_df = tpm_df[extreme_occurrences < max_extreme_occurence]
    tpm_df = tpm_df.fillna(0)
    print(f"✅ KO après nettoyage des valeurs extrêmes : {tpm_df.shape[0]}")


    #tpm_filtered = filter_tpm(tpm_df, abundance_threshold, presence_threshold, presence_fraction)
    tpm_filtered = filter_tpm(tpm_df, presence_fraction)
    print("💾 Sauvegarde du fichier TPM...")
    tpm_filtered.to_parquet(output_files[label])
    print(f"✅ {label.upper()} : {tpm_filtered.shape[0]} KO conservés après filtrage.")
