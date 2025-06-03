#!/usr/bin/env python3
"""
Build wide-format feature tables for mg and mt data from IMP output sample directories.

For each sample directory under <input_dir>, the script processes:
  - mg:
      - Taxonomic counts from Bracken files (located in Analysis/taxonomy/bracken/mg/).
      - Functional counts from annotation files matching pattern mg.*counts.tsv(.gz).
      - Library-size normalization is applied using mg-specific keys.
  - mt:
      - Functional counts from annotation files matching pattern mt.*counts.tsv(.gz).
      - Taxonomic features parsed from Kraken2 reports (located in Analysis/taxonomy/kraken/mt.kraken.report(.gz)).
      - For Kraken2, relative expression (column 1) is used.
      - Library-size normalization is applied using mt-specific keys.

Features are filtered based on a minimum value in at least a given fraction (prevalence) of samples.
Duplicate feature columns (with identical sanitized names) are summed.
Two Parquet tables (one for mg and one for mt) are saved using Snappy compression.

Usage Example:
    ./build_feature_tables.py --input_dir /path/to/imp/outs --out_dir /path/to/output

Dependencies:
    numpy, pandas, pyarrow, tqdm, scipy
"""

import os
import re
import glob
import gzip
import argparse
import logging
import concurrent.futures
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Tuple

# --- Constants and Precompiled Regex Patterns ---
_INVALID_CHARS_RE = re.compile(r"[^a-zA-Z0-9_]")
_MULTIPLE_UNDERSCORES_RE = re.compile(r"_+")
_LEADING_DIGIT_RE = re.compile(r"^\d")

# --- Logging Setup ---
class ExcludeWarningFilter(logging.Filter):
    """Filter to exclude WARNING level logs from a handler."""
    def filter(self, record: logging.LogRecord) -> bool:
        return record.levelno != logging.WARNING

def setup_logging(out_dir: str) -> None:
    """
    Configure logging for both console and file outputs.

    Args:
        out_dir (str): Directory where the log file will be saved.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.addFilter(ExcludeWarningFilter())
    console_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)

    # File handler (logs warnings and above)
    log_file = os.path.join(out_dir, "build_feature_tables.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.WARNING)
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

# --- Utility Functions ---
def smart_read_csv(path: str, **kwargs) -> pd.DataFrame:
    """
    Read a CSV file, supporting gzip compression.

    Args:
        path (str): File path.
        **kwargs: Additional keyword arguments for pandas.read_csv.

    Returns:
        pd.DataFrame: DataFrame containing the CSV data.
    """
    if path.endswith(".gz"):
        return pd.read_csv(path, compression="gzip", **kwargs)
    return pd.read_csv(path, **kwargs)

def sanitize_feature_name(name: str) -> str:
    """
    Sanitize a feature name by replacing invalid characters and ensuring it does not start with a digit.

    Args:
        name (str): Original feature name.

    Returns:
        str: Sanitized feature name.
    """
    name = _INVALID_CHARS_RE.sub("_", name)
    name = _MULTIPLE_UNDERSCORES_RE.sub("_", name)
    if _LEADING_DIGIT_RE.match(name):
        name = "f_" + name
    return name

def build_prefix(prefix: str, extra: str) -> str:
    """
    Build a feature name by prepending a prefix to a sanitized extra string.

    Args:
        prefix (str): The prefix string (e.g., 'mg' or 'mt').
        extra (str): Additional descriptor for the feature.

    Returns:
        str: Constructed feature name.
    """
    return f"{prefix}.{sanitize_feature_name(extra)}"

# --- Data Reading Functions ---
def read_functional_counts(func_path: str, prefix: str = "") -> tuple:
    """
    Load functional counts from an annotation file using NumPy.

    Reads column 0 (Geneid, string) and columns 6 and 7 (counts, float64), and returns
    the feature names and the sum of the counts.

    Args:
        func_path (str): Path to the annotation file.
        prefix (str): Prefix to add to feature names.

    Returns:
        tuple: (np.array of feature names, np.array of counts)
    """
    p = Path(func_path)
    if not p.exists():
        return np.array([], dtype="U50"), np.array([], dtype=np.float64)
    try:
        dt = np.dtype([("Geneid", "U50"), ("count1", np.float64), ("count2", np.float64)])
        data = np.genfromtxt(str(p), delimiter="\t", dtype=dt, comments="#", usecols=(0, 6, 7))
    except Exception as e:
        logging.error(f"Error reading {p}: {e}")
        return np.array([], dtype="U50"), np.array([], dtype=np.float64)
    if data.size == 0:
        logging.warning(f"No data found in {p}.")
        return np.array([], dtype="U50"), np.array([], dtype=np.float64)
    if data.ndim == 0:
        data = np.array([data])
    names = np.char.strip(data["Geneid"])
    vectorized_prefix = np.vectorize(lambda x: build_prefix(prefix, x))
    feats = vectorized_prefix(names)
    vals = data["count1"] + data["count2"]
    return feats, vals

def read_bracken(bracken_path: str, prefix: str = "") -> tuple:
    """
    Load taxonomic counts from a Bracken file using NumPy.

    Reads column 0 (feature name, string) and column 6 (fraction_total_reads, float32).

    Args:
        bracken_path (str): Path to the Bracken file.
        prefix (str): Prefix to add to feature names.

    Returns:
        tuple: (np.array of feature names, np.array of fractions)
    """
    p = Path(bracken_path)
    if not p.exists() and (p.with_suffix(p.suffix + ".gz")).exists():
        p = p.with_suffix(p.suffix + ".gz")
    if not p.exists():
        logging.error(f"Bracken file not found: {bracken_path}")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    try:
        dt = np.dtype([("name", "U50"), ("fraction", np.float32)])
        data = np.genfromtxt(str(p), delimiter="\t", dtype=dt,
                             comments="#", skip_header=1, usecols=(0, 6))
    except Exception as e:
        logging.error(f"Error reading Bracken file {p}: {e}")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    if data.size == 0:
        logging.warning(f"No data found in {p}.")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    if data.ndim == 0:
        data = np.array([data])
    names = np.char.strip(data["name"])
    vectorized_prefix = np.vectorize(lambda x: build_prefix(prefix, x))
    feats = vectorized_prefix(names)
    return feats, data["fraction"]

def read_kraken(kraken_path: str, prefix: str = "", tax_levels: list = None) -> tuple:
    """
    Parse a Kraken report file to extract taxonomic features and their relative expressions.

    Args:
        kraken_path (str): Path to the Kraken report.
        prefix (str): Prefix to add to feature names.
        tax_levels (list, optional): List of taxonomic levels to retain. Defaults to None.

    Returns:
        tuple: (np.array of feature names, np.array of relative expression values)
    """
    p = Path(kraken_path)
    if not p.exists() and (p.with_suffix(p.suffix + ".gz")).exists():
        p = p.with_suffix(p.suffix + ".gz")
    if not p.exists():
        logging.error(f"Kraken report not found: {kraken_path}")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    feats_list = []
    vals_list = []
    open_func = gzip.open if p.suffix == ".gz" else open
    try:
        with open_func(str(p), "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                tax_rank = parts[5].strip()
                if tax_levels is not None and tax_rank not in tax_levels:
                    continue
                try:
                    rel_expr = np.float32(parts[0])
                except ValueError:
                    continue
                taxon_name = parts[7].strip()
                feat = build_prefix(f"{prefix}.{tax_rank}", taxon_name)
                feats_list.append(feat)
                vals_list.append(rel_expr)
    except Exception as e:
        logging.error(f"Error reading Kraken report at {p}: {e}")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    if not feats_list:
        logging.warning(f"No data parsed from {p}.")
        return np.array([], dtype="U50"), np.array([], dtype=np.float32)
    return np.array(feats_list, dtype="U50"), np.array(vals_list, dtype=np.float32)

def get_library_size(stats_path: str, paired_key: str, paired_key2: str, single_key: str) -> int:
    """
    Compute the library size by summing read numbers from specified keys in a stats file.

    Args:
        stats_path (str): Path to the stats file.
        paired_key (str): Key for first paired read.
        paired_key2 (str): Key for second paired read.
        single_key (str): Key for single read.

    Returns:
        int: Total library size (read count).
    """
    p = Path(stats_path)
    if not p.exists():
        return 0
    try:
        open_func = gzip.open if p.suffix == ".gz" else open
        with open_func(str(p), "rt") as f:
            header = f.readline().strip().split("\t")
        try:
            step_idx = header.index("step")
            number_idx = header.index("number")
        except ValueError:
            logging.warning(f"Required columns 'step' or 'number' not found in {stats_path}.")
            return 0
        data = np.genfromtxt(str(p), delimiter="\t", dtype=None, encoding="utf-8",
                               names=True, usecols=(step_idx, number_idx))
    except Exception as e:
        logging.warning(f"Error reading stats file {stats_path}: {e}")
        return 0
    if data.size == 0:
        logging.warning(f"No data found in {stats_path}.")
        return 0
    if data.ndim == 0:
        data = np.array([data])
    keys = np.array([paired_key, paired_key2, single_key])
    mask = np.isin(data["step"], keys)
    if not np.any(mask):
        logging.warning(f"Could not find read numbers for {stats_path}.")
        return 0
    total = np.sum(data["number"][mask])
    return int(total)

# --- File-Finding Helper Functions ---
def find_bracken_file(sample_dir: str, ome: str, tax_level: str) -> str:
    """
    Recursively search for a Bracken file in the sample directory.

    Args:
        sample_dir (str): Path to the sample directory.
        ome (str): Omic type (e.g., 'mg').
        tax_level (str): Taxonomic level to search for.

    Returns:
        str: Path to the Bracken file, or an empty string if not found.
    """
    base_dir = Path(sample_dir) / "Analysis" / "taxonomy" / "bracken"
    patterns = [f"{ome}.{tax_level}*.bracken*", f"{tax_level}*.bracken*"]
    for pat in patterns:
        matches = list(base_dir.rglob(pat))
        if matches:
            if len(matches) > 1:
                logging.warning(f"Multiple Bracken files found for {ome} tax level {tax_level} using pattern {pat}: {mat
ches}. Using the first one.")
            return str(matches[0])
    try:
        contents = os.listdir(str(base_dir))
    except Exception as e:
        contents = f"Error reading directory: {e}"
    logging.error(f"No Bracken file found for {ome} tax level {tax_level} in {base_dir}. Directory contents: {contents}"
)
    return ""

def find_kraken_file(sample_dir: str, ome: str) -> str:
    """
    Recursively search for a Kraken report file in the sample directory.

    Args:
        sample_dir (str): Path to the sample directory.
        ome (str): Omic type (e.g., 'mt').

    Returns:
        str: Path to the Kraken report file, or an empty string if not found.
    """
    base_dir = Path(sample_dir) / "Analysis" / "taxonomy" / "kraken"
    patterns = [f"{ome}.kraken.report*", "kraken.report*"]
    for pat in patterns:
        matches = list(base_dir.rglob(pat))
        if matches:
            if len(matches) > 1:
                logging.warning(f"Multiple Kraken report files found for {ome} using pattern {pat}: {matches}. Using the
 first one.")
            return str(matches[0])
    try:
        contents = os.listdir(str(base_dir))
    except Exception as e:
        contents = f"Error reading directory: {e}"
    logging.error(f"No Kraken report file found for {ome} in {base_dir}. Directory contents: {contents}")
    return ""

# --- Data Gathering and Processing Functions ---
def gather_sample_data(sample_dir: str, ome: str, tax_levels: list, stats_filename: str,
                       paired_key: str, paired_key2: str, single_key: str) -> dict:
    """
    Read and aggregate taxonomic and functional counts for a given sample.

    Args:
        sample_dir (str): Path to the sample directory.
        ome (str): Omic type ('mg' or 'mt').
        tax_levels (list): List of taxonomic levels to process.
        stats_filename (str): Stats file name.
        paired_key (str): Paired read key.
        paired_key2 (str): Second paired read key.
        single_key (str): Single read key.

    Returns:
        dict: Dictionary mapping sanitized feature names to normalized counts.
    """
    sample = os.path.basename(sample_dir)
    stats_path = os.path.join(sample_dir, "Stats", stats_filename)
    lib_size = get_library_size(stats_path, paired_key, paired_key2, single_key)
    out_data = {}

    # Process taxonomic data
    if ome == "mg":
        for lvl in tqdm(tax_levels, desc=f"Processing {sample} {ome} tax levels", leave=False, unit="level", position=1)
:
            bracken_path = find_bracken_file(sample_dir, ome, lvl)
            if bracken_path:
                prefix = f"{ome}.{lvl}"
                feats, vals = read_bracken(bracken_path, prefix)
                if feats.size > 0:
                    for feat, val in zip(feats, vals):
                        out_data[feat] = val
    elif ome == "mt":
        kraken_path = find_kraken_file(sample_dir, ome)
        if kraken_path:
            feats, vals = read_kraken(kraken_path, prefix=ome, tax_levels=tax_levels)
            for feat, val in zip(feats, vals):
                out_data[feat] = val
        else:
            logging.warning(f"Kraken report not found for sample {sample_dir}.")

    # Process functional counts (annotations)
    annotation_dir = os.path.join(sample_dir, "Analysis", "annotation")
    pattern = os.path.join(annotation_dir, f"{ome}.*counts.tsv*")
    func_files = [f for f in glob.glob(pattern)
                  if not f.endswith("summary.gz") and not f.endswith("summary")]

    annotation_bar = tqdm(func_files, desc=f"Processing {sample} {ome} annotation files", leave=False, unit="file", posi
tion=1)
    for func_path in annotation_bar:
        if os.path.getsize(func_path) == 0:
            logging.warning(f"Annotation file {func_path} is empty. Skipping.")
            continue
        base = os.path.basename(func_path)
        if base.endswith(".gz"):
            base = base[:-3]
        if base.endswith(".tsv"):
            base = base[:-4]
        if base.startswith(f"{ome}."):
            base = base[len(ome) + 1:]
        prefix = f"{ome}.{base}"
        feats, vals = read_functional_counts(func_path, prefix)
        if lib_size > 0 and vals.size > 0:
            vals = vals / lib_size
        elif lib_size == 0:
            logging.warning(f"Library size is zero for sample {sample_dir}. Skipping normalization for {func_path}.")
        for feat, val in zip(feats, vals):
            out_data[feat] = val
    annotation_bar.close()

    # Sanitize feature names
    sanitized_data = {sanitize_feature_name(key): value for key, value in out_data.items()}
    return sanitized_data

def process_sample(sample: str, input_dir: str, ome: str, tax_levels: list, stats_filename: str,
                   paired_key: str, paired_key2: str, single_key: str) -> tuple:
    """
    Process a single sample directory and extract feature data.

    Args:
        sample (str): Sample identifier.
        input_dir (str): Base directory containing sample subdirectories.
        ome (str): Omic type ('mg' or 'mt').
        tax_levels (list): List of taxonomic levels.
        stats_filename (str): Stats file name.
        paired_key (str): Paired read key.
        paired_key2 (str): Second paired read key.
        single_key (str): Single read key.

    Returns:
        tuple: (sample ID, np.array of feature names, np.array of corresponding values)
    """
    sample_dir = os.path.join(input_dir, sample)
    features_dict = gather_sample_data(sample_dir, ome, tax_levels, stats_filename,
                                       paired_key, paired_key2, single_key)
    if not features_dict:
        return sample, np.array([], dtype="U50"), np.array([], dtype=np.float64)
    feat_arr = np.array(list(features_dict.keys()), dtype="U50")
    val_arr = np.array(list(features_dict.values()), dtype=np.float64)
    return sample, feat_arr, val_arr

def sort_columns(columns: list, ome: str, tax_levels: list) -> list:
    """
    Sort feature table columns so that taxonomic features appear first followed by annotation features.

    Args:
        columns (list): List of feature column names.
        ome (str): Omic type ('mg' or 'mt').
        tax_levels (list): Order of taxonomic levels.

    Returns:
        list: Sorted list of columns.
    """
    def col_key(col: str) -> tuple:
        parts = col.split(".")
        if len(parts) >= 3:
            if parts[1] in tax_levels:
                tax_order = tax_levels.index(parts[1])
                tax_name = parts[2]
                remainder = ".".join(parts[3:]) if len(parts) > 3 else ""
                return (0, tax_order, tax_name, remainder, col)
            else:
                block = parts[1]
                remainder = ".".join(parts[2:])
                return (1, block, remainder, col)
        return (2, col)
    return sorted(columns, key=col_key)

def build_feature_table(input_dir: str, ome: str, tax_levels: list, stats_filename: str,
                        paired_key: str, paired_key2: str, single_key: str,
                        normalization_threshold: float, prevalence: float,
                        num_workers: int, clr_transform: bool = False) -> pd.DataFrame:
    """
    Build a wide-format feature table for the given omic type.

    The table is constructed using dense NumPy arrays and filtered based on feature prevalence.
    Optionally applies centered log-ratio (CLR) transformation.

    Args:
        input_dir (str): Base input directory.
        ome (str): Omic type ('mg' or 'mt').
        tax_levels (list): Taxonomic levels to process.
        stats_filename (str): Stats file name.
        paired_key (str): Paired read key.
        paired_key2 (str): Second paired read key.
        single_key (str): Single read key.
        normalization_threshold (float): Minimum value threshold for a feature.
        prevalence (float): Minimum fraction of samples in which a feature must be present.
        num_workers (int): Number of worker processes (-1 uses all available cores).
        clr_transform (bool): If True, apply CLR transformation.

    Returns:
        pd.DataFrame: DataFrame containing the feature table.
    """
    samples = [d for d in os.listdir(input_dir)
               if os.path.isdir(os.path.join(input_dir, d)) and not d.startswith(".")]
    logging.info(f"Found {len(samples)} sample directories for {ome} in {input_dir}.")

    sample_data = []
    if num_workers == -1:
        num_workers = os.cpu_count() or 1

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(process_sample, sample, input_dir, ome, tax_levels, stats_filename,
                            paired_key, paired_key2, single_key): sample for sample in samples
        }
        with tqdm(total=len(futures), desc=f"Processing {ome} samples", unit="sample") as pbar:
            for future in concurrent.futures.as_completed(futures):
                try:
                    sample_data.append(future.result())
                except Exception as e:
                    logging.error(f"Error processing sample {futures[future]}: {e}")
                pbar.update(1)
    if not sample_data:
        logging.error(f"No {ome} sample records were processed.")
        return pd.DataFrame()

    sample_ids = [sample for sample, feats, vals in sample_data]
    logging.info(f"Processed {len(sample_ids)} samples.")

    # Build union of features
    logging.info("Building union of all features...")
    all_feats = np.concatenate([feats for _, feats, _ in sample_data if feats.size > 0])
    unique_features = np.unique(all_feats)
    logging.info(f"Found {len(unique_features)} unique features.")

    n_samples = len(sample_ids)
    n_features = unique_features.shape[0]
    out_array = np.zeros((n_samples, n_features), dtype=np.float32)

    # Populate dense matrix
    for i, (sample, feats, vals) in enumerate(tqdm(sample_data, desc="Filling dense matrix", unit="sample")):
        if feats.size == 0:
            continue
        col_indices = np.searchsorted(unique_features, feats)
        np.add.at(out_array[i], col_indices, vals)

    # Prevalence filtering
    col_counts = (out_array > normalization_threshold).sum(axis=0)
    min_count = int(prevalence * n_samples)
    logging.info(f"Filtering: retaining features present in at least {min_count} samples ({prevalence*100:.1f}%).")
    keep = col_counts >= min_count
    filtered_features = unique_features[keep]
    out_array = out_array[:, keep]
    logging.info(f"Retained {out_array.shape[1]} features after filtering.")

    # Optional CLR transformation
    if clr_transform:
        arr_pos = np.where(out_array > 0, out_array, np.inf)
        min_pos = np.min(arr_pos, axis=0)
        pseudocount = min_pos * 0.65
        gm = np.exp(np.mean(np.log(out_array + pseudocount), axis=1))
        out_array = np.log((out_array + pseudocount) / gm[:, None]).astype(np.float32)

    df = pd.DataFrame(out_array, index=sample_ids, columns=filtered_features)
    sorted_cols = sort_columns(list(df.columns), ome, tax_levels)
    df = df[sorted_cols]

    return df

def check_ome_exists(input_dir: str, ome: str) -> bool:
    """
    Check if preprocessing files exist for the specified omic type in any sample.

    Args:
        input_dir (str): Base input directory.
        ome (str): Omic type ('mg' or 'mt').

    Returns:
        bool: True if at least one sample contains preprocessing files for the omic type.
    """
    samples = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    for sample in samples:
        sample_dir = os.path.join(input_dir, sample)
        preprocessing_dir = os.path.join(sample_dir, "Preprocessing")
        if os.path.isdir(preprocessing_dir):
            pattern = os.path.join(preprocessing_dir, f"{ome}*")
            if glob.glob(pattern):
                return True
    return False

def write_feature_summary(df: pd.DataFrame, summary_file: str, tax_levels: list, ome: str) -> None:
    """
    Write a long-format TSV summary file listing each feature along with its feature group and source.

    Taxonomic features (those without '_counts') are labeled as "taxonomy" and are ordered by the provided
    tax_levels (from highest to lowest rank). Annotation features (those with '_counts') are labeled as "annotation"
    and sorted alphabetically by feature group and then by feature name.

    The overall sorting order is defined as follows:
        1. Taxonomic features precede annotation features.
        2. Taxonomic features are sorted according to tax_levels (high levels first).
        3. Within each group, features are sorted alphabetically by the feature group.
        4. Finally, features are sorted alphabetically by their full name.

    Args:
        df (pd.DataFrame): DataFrame containing the feature table.
        summary_file (str): Output path for the summary file.
        tax_levels (list): List of taxonomic levels in order from highest to lowest (e.g., ["R1", "D", "P", "C", "O", "F
", "G", "S"]).
        ome (str): Omic type ('mg' or 'mt').
    """

    def extract_feature_info(col: str, ome: str) -> Tuple[str, str]:
        """
        Extract the feature group and source from a column name.

        Taxonomic features are expected to not include '_counts' and begin with the omic prefix.
        Annotation features include '_counts'.

        Args:
            col (str): The column name.
            ome (str): Omic type.

        Returns:
            Tuple[str, str]: A tuple of (feature_group, source) where source is either "taxonomy" or "annotation".
        """
        col_us = col.replace(".", "_")
        prefix = f"{ome}_"
        remainder = col_us[len(prefix):] if col_us.startswith(prefix) else col_us
        if "_counts" in remainder:
            group = remainder.split("_counts")[0]
            source = "annotation"
        else:
            # For taxonomy, assume the group is the prefix up to the first underscore.
            parts = remainder.split("_", 1)
            group = parts[0] if parts else remainder
            source = "taxonomy"
        return group, source

    # Build a mapping from taxonomic level to its order (lower values indicate higher rank).
    taxonomy_order = {lvl: i for i, lvl in enumerate(tax_levels)}

    rows = []
    for col in df.columns:
        group, source = extract_feature_info(col, ome)
        if source == "taxonomy":
            # Get order based on taxonomy_order; default to a high number if group is not found.
            level_order = taxonomy_order.get(group, len(tax_levels))
            sort_key = (0, level_order, group.lower(), col.lower())
        else:
            # Annotation features follow taxonomy features.
            sort_key = (1, 0, group.lower(), col.lower())
        rows.append({
            "feature": col,
            "feature_group": group,
            "source": source,
            "sort_key": sort_key
        })

    df_summary = pd.DataFrame(rows)
    df_summary.sort_values(by="sort_key", inplace=True)
    df_summary.drop(columns=["sort_key"], inplace=True)
    df_summary.to_csv(summary_file, sep="\t", index=False)

# --- Main Execution Function ---
def main() -> None:
    """
    Main function to parse command-line arguments and build mg and mt feature tables.
    """
    parser = argparse.ArgumentParser(
        description="Build separate wide-format feature tables for mg and mt data from IMP sample directories."
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing sample subdirectories (e.g., path/to/imp/outs).")
    parser.add_argument("--out_dir", required=True,
                        help="Output directory to store the mg and mt feature table Parquet files.")
    parser.add_argument("--normalization_threshold", type=float, default=1e-7,
                        help="Minimum value a feature must have in a sample (default: 1e-7).")
    parser.add_argument("--prevalence", type=float, default=0.3,
                        help="Minimum fraction of samples that must exceed the threshold (default: 0.3).")
    parser.add_argument("--tax_levels", nargs="*", default=["D", "P", "G", "S"],
                        help="Taxonomic levels to process (default: D P G S).")
    parser.add_argument("--stats_filename", default="all_stats.tsv.gz",
                        help="Stats filename (default: all_stats.tsv.gz).")
    parser.add_argument("--num_workers", type=int, default=-1,
                        help="Number of worker processes to use (-1 means use all available cores).")
    # mg-specific keys.
    parser.add_argument("--mg_paired_key", default="mg.r1.trimmed__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v2
transcripts_filtered.fq",
                        help="Paired read key for mg.")
    parser.add_argument("--mg_paired_key2", default="mg.r2.trimmed__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v
2transcripts_filtered.fq",
                        help="Second paired read key for mg.")
    parser.add_argument("--mg_single_key", default="mg.se.trimmed__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v2
transcripts_filtered.fq",
                        help="Single read key for mg.")
    # mt-specific keys.
    parser.add_argument("--mt_paired_key", default="mt.r1.trimmed__rrna_proc__PhiX_filtered__T2TCHM13v2genomic_filtered_
_T2TCHM13v2transcripts_filtered.fq",
                        help="Paired read key for mt.")
    parser.add_argument("--mt_paired_key2", default="mt.r2.trimmed__rrna_proc__PhiX_filtered__T2TCHM13v2genomic_filtered
__T2TCHM13v2transcripts_filtered.fq",
                        help="Second paired read key for mt.")
    parser.add_argument("--mt_single_key", default="mt.se.trimmed__rrna_proc__PhiX_filtered__T2TCHM13v2genomic_filtered_
_T2TCHM13v2transcripts_filtered.fq",
                        help="Single read key for mt.")
    parser.add_argument("--clr_transform", action="store_true", default=False,
                        help="If set, apply centered log-ratio (CLR) transformation to the feature table.")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    setup_logging(args.out_dir)
    logging.info("Building mg and mt feature tables from sample directories...")

    # Process mg features
    if check_ome_exists(args.input_dir, "mg"):
        df_mg = build_feature_table(
            input_dir=args.input_dir,
            ome="mg",
            tax_levels=args.tax_levels,
            stats_filename=args.stats_filename,
            paired_key=args.mg_paired_key,
            paired_key2=args.mg_paired_key2,
            single_key=args.mg_single_key,
            normalization_threshold=args.normalization_threshold,
            prevalence=args.prevalence,
            num_workers=args.num_workers,
            clr_transform=args.clr_transform
        )
        if df_mg.empty:
            logging.error("No mg features were processed. Exiting mg processing.")
        else:
            mg_out_path = os.path.join(args.out_dir, "mg_feature_table.parquet")
            df_mg = df_mg.reset_index()  # "index" holds sample IDs.
            df_mg.to_parquet(mg_out_path, engine="pyarrow", compression="snappy")
            logging.info(f"mg feature table saved to {mg_out_path}")
            mg_summary_path = os.path.join(args.out_dir, "mg_feature_summary.tsv")
            feature_df = df_mg.drop(columns=["index"])
            write_feature_summary(feature_df, mg_summary_path, tax_levels=args.tax_levels, ome="mg")
            logging.info(f"mg feature summary saved to {mg_summary_path}")
    else:
        logging.error("No mg Preprocessing files found in any sample. Skipping mg feature table construction.")

    # Process mt features
    if check_ome_exists(args.input_dir, "mt"):
        df_mt = build_feature_table(
            input_dir=args.input_dir,
            ome="mt",
            tax_levels=args.tax_levels,
            stats_filename=args.stats_filename,
            paired_key=args.mt_paired_key,
            paired_key2=args.mt_paired_key2,
            single_key=args.mt_single_key,
            normalization_threshold=args.normalization_threshold,
            prevalence=args.prevalence,
            num_workers=args.num_workers,
            clr_transform=args.clr_transform
        )
        if df_mt.empty:
            logging.error("No mt features were processed. Exiting mt processing.")
        else:
            mt_out_path = os.path.join(args.out_dir, "mt_feature_table.parquet")
            df_mt = df_mt.reset_index()  # "index" holds sample IDs.
            df_mt.to_parquet(mt_out_path, engine="pyarrow", compression="snappy")
            logging.info(f"mt feature table saved to {mt_out_path}")
            mt_summary_path = os.path.join(args.out_dir, "mt_feature_summary.tsv")
            feature_df = df_mt.drop(columns=["index"])
            write_feature_summary(feature_df, mt_summary_path, tax_levels=args.tax_levels, ome="mt")
            logging.info(f"mt feature summary saved to {mt_summary_path}")
    else:
        logging.error("No mt Preprocessing files found in any sample. Skipping mt feature table construction.")

if __name__ == "__main__":
    main()
