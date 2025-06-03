import pyarrow.parquet as pq
import pyarrow as pa

from tabulate import tabulate


file_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/mt_feature_table.parquet"

pq_array = pq.read_table(file_path, memory_map=True)
print("RSS with memory_map=True: {}MB".format(pa.total_allocated_bytes() >> 20))

pq_array = pq.read_table(file_path, memory_map=False)
print("RSS with memory_map=False: {}MB".format(pa.total_allocated_bytes() >> 20))

print(pq_array.to_pandas().head())

print(pq_array.to_pandas().head(15))  # Afficher les 5 premières lignes
