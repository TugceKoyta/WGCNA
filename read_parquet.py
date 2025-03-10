import pyarrow.parquet as pq
import pyarrow as pa


# Chemin du fichier Parquet (à modifier)
file_path = "/scratch/users/tkoytaviloglu/results/outputs/mt_feature_table.parquet"

# Lire le fichier avec memory_map=True
pq_array = pq.read_table(file_path, memory_map=True)
print("RSS with memory_map=True: {}MB".format(pa.total_allocated_bytes() >> 20))

# Lire le fichier avec memory_map=False
pq_array = pq.read_table(file_path, memory_map=False)
print("RSS with memory_map=False: {}MB".format(pa.total_allocated_bytes() >> 20))

# Afficher un aperçu des données
print(pq_array.to_pandas().head())

print(pq_array.to_pandas().head())  # Afficher les 5 premières lignes
