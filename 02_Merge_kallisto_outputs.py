"""
Merge Kallisto 'abundance.tsv' TPM columns across samples into a single table.

Features:
- Uses argparse for runtime configuration
- Robust checking of per-sample abundance files
- Efficient concatenation using pandas.concat on Series
- Optional transcript->gene mapping join (GTF/TSV) if provided
- Options to fill missing values, compress output, and keep logs

Example usage:
python cleaned_kallisto_merge.py \
  --kallisto-dir /path/to/kallisto_fastq_output/tumor/ \
  --output merged_tpm.tsv.gz \
  --abundance-file-name abundance.tsv \
  --fill-na 0

"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import logging
import sys
from typing import Optional


def setup_logger(log_file: Optional[Path] = None):
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[logging.StreamHandler(sys.stdout)])
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter(fmt))
        logging.getLogger().addHandler(fh)


def find_sample_dirs(base_dir: Path, depth_one: bool = True):
    """Return subdirectories directly under base_dir (one level by default).
    If depth_one is False, search recursively for directories that contain an abundance file.
    """
    if depth_one:
        return sorted([p for p in base_dir.iterdir() if p.is_dir()])
    else:
        # recursive search: directories that contain abundance files
        return sorted({p.parent for p in base_dir.rglob("*") if p.is_dir() or p.name})


def read_abundance_tpm(abundance_path: Path, id_col: str = "target_id", tpm_col: str = "tpm") -> Optional[pd.Series]:
    """Read an abundance.tsv and return a pandas Series indexed by id_col with values from tpm_col.
    Returns None if file missing or malformed.
    """
    if not abundance_path.exists():
        logging.warning("Missing file: %s", abundance_path)
        return None
    try:
        df = pd.read_csv(abundance_path, sep="\t", usecols=[id_col, tpm_col], dtype={id_col: str})
    except ValueError as e:
        logging.warning("File %s missing required columns (%s). Skipping.", abundance_path, e)
        return None
    except pd.errors.EmptyDataError:
        logging.warning("File %s is empty. Skipping.", abundance_path)
        return None

    if df.empty:
        logging.warning("File %s contains no rows. Skipping.", abundance_path)
        return None

    # In case of duplicated target_id, aggregate (sum) TPMs for same id (unlikely but safe)
    if df[id_col].duplicated().any():
        logging.debug("Found duplicated %s in %s; aggregating by sum.", id_col, abundance_path)
        df = df.groupby(id_col, as_index=False)[tpm_col].sum()

    s = pd.Series(df[tpm_col].values, index=df[id_col].values, name=abundance_path.parent.name)
    s.index.name = id_col
    return s


def merge_kallisto_tpm(
    base_dir: Path,
    abundance_file_name: str = "abundance.tsv",
    id_col: str = "target_id",
    tpm_col: str = "tpm",
    recursive: bool = False,
):
    """Scan sample subdirectories and merge TPM columns into a DataFrame.

    Returns the merged DataFrame (index = id_col, columns = sample names).
    """
    base_dir = Path(base_dir)
    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory does not exist: {base_dir}")

    # Get candidate sample directories
    if recursive:
        # find directories containing abundance files
        sample_dirs = sorted([p.parent for p in base_dir.rglob(abundance_file_name) if p.is_file()])
    else:
        sample_dirs = sorted([p for p in base_dir.iterdir() if p.is_dir()])

    logging.info("Found %d candidate sample directories", len(sample_dirs))

    series_list = []
    sample_names = []

    for d in sample_dirs:
        abundance_path = d / abundance_file_name
        s = read_abundance_tpm(abundance_path, id_col=id_col, tpm_col=tpm_col)
        if s is None:
            continue
        series_list.append(s)
        sample_names.append(s.name)

    if len(series_list) == 0:
        raise RuntimeError("No abundance files loaded. Check the directory and file name pattern.")

    # Efficient concatenation: align on the index (target_id)
    merged = pd.concat(series_list, axis=1, sort=True)
    merged.index.name = id_col
    return merged


def optionally_add_gene_map(merged_df: pd.DataFrame, mapping_file: Optional[Path], id_col: str = "target_id") -> pd.DataFrame:
    """If mapping_file is provided (TSV with columns: target_id and gene_id or gene_name), join gene info.
    mapping_file should have at least two columns, the first matching target_id and the second with gene names/ids.
    """
    if mapping_file is None:
        return merged_df
    mapping_file = Path(mapping_file)
    if not mapping_file.exists():
        logging.warning("Mapping file not found: %s. Skipping gene mapping.", mapping_file)
        return merged_df

    try:
        map_df = pd.read_csv(mapping_file, sep="\t", dtype={id_col: str}, header=0)
    except Exception as e:
        logging.warning("Could not read mapping file %s: %s. Skipping mapping.", mapping_file, e)
        return merged_df

    # Guess which columns to use
    cols = map_df.columns.tolist()
    if len(cols) < 2:
        logging.warning("Mapping file has fewer than 2 columns; expected at least target_id and gene name. Skipping.")
        return merged_df

    map_df = map_df.rename(columns={cols[0]: id_col})
    map_df = map_df.drop_duplicates(subset=id_col)
    map_df.set_index(id_col, inplace=True)

    # Left-join mapping to merged_df (mapping columns appended on the right)
    out = merged_df.join(map_df, how="left")
    return out


def parse_args():
    p = argparse.ArgumentParser(description="Merge Kallisto abundance TPMs across samples")
    p.add_argument("--kallisto-dir", required=True, help="Base directory containing sample subdirectories with abundance.tsv")
    p.add_argument("--abundance-file-name", default="abundance.tsv", help="Name of the abundance file inside each sample dir")
    p.add_argument("--output", required=True, help="Path to write merged output (tsv or tsv.gz supported)")
    p.add_argument("--id-col", default="target_id", help="Column name for transcript/gene id in abundance file")
    p.add_argument("--tpm-col", default="tpm", help="Column name for TPM values in abundance file")
    p.add_argument("--recursive", action="store_true", help="Search recursively for abundance files rather than only direct subdirectories")
    p.add_argument("--map-file", default=None, help="Optional TSV mapping file (target_id -> gene_name) to append columns")
    p.add_argument("--fill-na", type=str, default=None, help="Value to fill missing TPMs (e.g. 0). Use 'None' to keep NaN")
    p.add_argument("--log-file", default=None, help="Optional path to a log file")
    return p.parse_args()


def main():
    args = parse_args()
    setup_logger(args.log_file)

    kallisto_dir = Path(args.kallisto_dir)
    out_path = Path(args.output)

    merged = merge_kallisto_tpm(
        kallisto_dir,
        abundance_file_name=args.abundance_file_name,
        id_col=args.id_col,
        tpm_col=args.tpm_col,
        recursive=args.recursive,
    )

    # Optionally add mapping
    merged = optionally_add_gene_map(merged, Path(args.map_file) if args.map_file else None, id_col=args.id_col)

    # Fill NA if requested
    if args.fill_na is not None and args.fill_na.lower() != "none":
        try:
            fill_value = float(args.fill_na)
        except ValueError:
            fill_value = args.fill_na
        merged = merged.fillna(fill_value)

    # Save output (automatically compress if .gz extension provided)
    if out_path.suffix == ".gz":
        merged.to_csv(out_path, sep="\t", compression="gzip")
    else:
        merged.to_csv(out_path, sep="\t")

    logging.info("Merged data saved to %s", out_path)


if __name__ == "__main__":
    main()
