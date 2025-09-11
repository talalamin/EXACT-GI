"""
GTF extractor for a Imprinted gene list

This script extracts gene and/or transcript records from a GTF (or GTF.gz) file
for genes listed in a plain text file. It writes a compact TSV with selected

- CLI (argparse) for inputs and options
- Robust GTF attribute parsing
- Supports gzipped GTF files transparently
- Options to include gene and/or transcript features
- Optionally include additional attribute columns if present

Example:
python 03_gtf_gene_extractor.py \
  --gene-list tine_imprinted_genes.list \
  --gtf Homo_sapiens.GRCh38.108.gtf.gz \
  --out GI_tine_genes_info.tsv.gz \
  --features gene transcript

"""

from __future__ import annotations
import argparse
from pathlib import Path
import gzip
import io
import logging
import sys
from typing import Dict, Iterable, List, Optional


def setup_logger(log_file: Optional[Path] = None):
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[logging.StreamHandler(sys.stdout)])
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter(fmt))
        logging.getLogger().addHandler(fh)


def read_gene_list(path: Path) -> List[str]:
    """Read gene names from a text file. Strips quotes and whitespace.
    One gene name per line is expected. Empty lines are ignored.
    """
    gene_names = []
    with path.open("r", encoding="utf8") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            # remove surrounding quotes if present
            if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
                s = s[1:-1]
            gene_names.append(s)
    # unique and preserve order
    seen = set()
    out = []
    for g in gene_names:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def open_text(path: Path):
    """Open plain or gzipped text file for reading, return a text-mode file object."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, mode="rb"), encoding="utf8")
    return path.open("r", encoding="utf8")


def parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    """Parse the attribute column from a GTF/GFF3-like line into a dict.

    GTF attributes are semi-colon separated key value pairs like:
    gene_id "ENSG..."; gene_name "TP53"; transcript_id "ENST...";

    This function is robust to irregular spacing and missing trailing semicolons.
    """
    attrs: Dict[str, str] = {}
    if not attr_string:
        return attrs
    # Split on semicolons first
    parts = [p.strip() for p in attr_string.strip().split(';') if p.strip()]
    for part in parts:
        # GTF usually: key "value" or key "value"; but be forgiving
        if ' ' in part:
            key, val = part.split(' ', 1)
            val = val.strip().strip('"')
            attrs[key] = val
        else:
            # fallback: if no space, try split on '=' (some GFFs use key=value)
            if '=' in part:
                key, val = part.split('=', 1)
                attrs[key.strip()] = val.strip().strip('"')
            else:
                # unknown format, store raw
                attrs[part] = ''
    return attrs


def iterate_gtf_matches(gtf_path: Path, gene_set: Iterable[str], features: List[str]):
    """Yield parsed records from GTF whose gene_name matches entries in gene_set.

    Yields a tuple of (gene_name, record_dict) where record_dict contains:
    chrom, source, feature, start, end, score, strand, frame, attributes(dict)
    """
    gene_lookup = set(gene_set)
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes_field = parts[:9]
            if feature not in features:
                continue
            attr = parse_gtf_attributes(attributes_field)
            # gene_name is common but might be missing; try gene_id if needed
            gene_name = attr.get('gene_name') or attr.get('gene_id')
            if gene_name and gene_name in gene_lookup:
                rec = {
                    'chrom': chrom,
                    'source': source,
                    'feature': feature,
                    'start': start,
                    'end': end,
                    'score': score,
                    'strand': strand,
                    'frame': frame,
                    'attributes': attr,
                }
                yield gene_name, rec


def write_results_tsv(out_path: Path, rows: List[List[str]], header: List[str], compress: bool = False):
    mode = 'wt'
    if compress or str(out_path).endswith('.gz'):
        with gzip.open(out_path, 'wt', encoding='utf8') as fh:
            fh.write('\t'.join(header) + '\n')
            for r in rows:
                fh.write('\t'.join(r) + '\n')
    else:
        with out_path.open('w', encoding='utf8') as fh:
            fh.write('\t'.join(header) + '\n')
            for r in rows:
                fh.write('\t'.join(r) + '\n')


def main():
    p = argparse.ArgumentParser(description='Extract gene/transcript info from GTF for a provided gene list')
    p.add_argument('--gene-list', required=True, help='Plain text file with one gene name per line')
    p.add_argument('--gtf', required=True, help='Path to GTF file (supports .gz)')
    p.add_argument('--out', required=True, help='Output TSV path (.tsv or .tsv.gz)')
    p.add_argument('--features', nargs='+', default=['gene', 'transcript'], help='Which features to include (e.g. gene transcript exon)')
    p.add_argument('--extra-attr', nargs='*', default=[], help='Additional attribute keys to include as columns (e.g. gene_version)')
    p.add_argument('--log-file', default=None, help='Optional log file path')
    args = p.parse_args()

    setup_logger(Path(args.log_file) if args.log_file else None)

    gene_list_path = Path(args.gene_list)
    gtf_path = Path(args.gtf)
    out_path = Path(args.out)

    if not gene_list_path.exists():
        logging.error('Gene list file not found: %s', gene_list_path)
        sys.exit(1)
    if not gtf_path.exists():
        logging.error('GTF file not found: %s', gtf_path)
        sys.exit(1)

    genes = read_gene_list(gene_list_path)
    logging.info('Loaded %d gene names from %s', len(genes), gene_list_path)

    features = args.features
    extra_attrs = args.extra_attr

    rows = []
    seen = set()

    # Header columns: Gene Name, Gene ID, Transcript ID, [extra attrs...], Strand, Coordinates, Feature, Source
    header = ['Gene Name', 'Gene ID', 'Transcript ID'] + [a for a in extra_attrs] + ['Strand', 'Coordinates', 'Feature', 'Source']

    for gene_name, rec in iterate_gtf_matches(gtf_path, genes, features):
        attr = rec['attributes']
        gene_id = attr.get('gene_id', '')
        transcript_id = attr.get('transcript_id', '')
        extra_vals = [attr.get(a, '') for a in extra_attrs]
        coords = f"{rec['chrom']}:{rec['start']}-{rec['end']}"
        row = [gene_name, gene_id, transcript_id] + extra_vals + [rec['strand'], coords, rec['feature'], rec['source']]
        # Avoid exact duplicate rows
        row_key = '\t'.join(row)
        if row_key in seen:
            continue
        seen.add(row_key)
        rows.append(row)

    if len(rows) == 0:
        logging.warning('No matching records found for provided gene list.')

    compress_out = str(out_path).endswith('.gz')
    write_results_tsv(out_path, rows, header, compress=compress_out)
    logging.info('Wrote %d rows to %s', len(rows), out_path)


if __name__ == '__main__':
    main()
