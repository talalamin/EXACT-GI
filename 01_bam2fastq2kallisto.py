"""
BAM -> FASTQ -> Kallisto pipeline
Cleaned, parameterized and documented script suitable for a public repository.

Features:
- Uses argparse for runtime configuration
- Supports reading a BAM list file or scanning a directory
- Option to use GATK SamToFastq or samtools fastq
- Parallel processing with a configurable worker pool
- Robust logging and error handling
- Option to keep or delete intermediate FASTQ files

Example usage:
python cleaned_bam_to_kallisto.py \
  --bam-list /path/to/bam_list.txt \
  --fastq-out /path/to/fastq_out \
  --kallisto-out /path/to/kallisto_out \
  --kallisto-index /path/to/kallisto_index.idx \
  --num-processes 8 \
  --kallisto-threads 4

"""

from __future__ import annotations
import argparse
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import logging
import sys

# -------------------------- Utility helpers -------------------------------

def setup_logger(log_file: Path | str | None = None):
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    logging.basicConfig(level=logging.INFO, format=fmt, handlers=[logging.StreamHandler(sys.stdout)])
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter(fmt))
        logging.getLogger().addHandler(fh)


def which_or_path(name_or_path: str) -> str | None:
    """Return executable path if available either by absolute path or PATH lookup."""
    if os.path.isabs(name_or_path) and os.access(name_or_path, os.X_OK):
        return name_or_path
    return shutil.which(name_or_path)


# -------------------------- Core processing -------------------------------

def run_cmd(cmd: list[str], dry_run: bool = False) -> tuple[int, str, str]:
    """Run a command and return (returncode, stdout, stderr)."""
    logging.debug("Running command: %s", " ".join(cmd))
    if dry_run:
        logging.info("Dry-run: %s", " ".join(cmd))
        return 0, "", ""
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def bam_to_fastq_and_kallisto(
    bam_path: str,
    fastq_out_dir: str,
    kallisto_out_dir: str,
    kallisto_index: str,
    use_gatk: bool = True,
    gatk_exec: str | None = "gatk",
    samtools_exec: str | None = "samtools",
    kallisto_exec: str | None = "kallisto",
    kallisto_threads: int = 4,
    keep_fastq: bool = False,
    dry_run: bool = False,
) -> dict:
    """
    Process a single BAM: convert to paired FASTQ and run kallisto quant.
    Returns a dict with status and messages.
    """
    bam_path = str(bam_path)
    sample_name = Path(bam_path).stem
    fastq_1 = os.path.join(fastq_out_dir, f"{sample_name}_1.fastq")
    fastq_2 = os.path.join(fastq_out_dir, f"{sample_name}_2.fastq")
    kallisto_sample_out = os.path.join(kallisto_out_dir, sample_name)

    os.makedirs(fastq_out_dir, exist_ok=True)
    os.makedirs(kallisto_sample_out, exist_ok=True)

    # 1) Convert BAM -> FASTQ
    if use_gatk:
        gatk_bin = which_or_path(gatk_exec or "gatk")
        if not gatk_bin:
            return {"sample": sample_name, "ok": False, "message": "GATK not found"}
        cmd_gatk = [gatk_bin, "SamToFastq", "-I", bam_path, "-F", fastq_1, "--SECOND_END_FASTQ", fastq_2]
        rc, out, err = run_cmd(cmd_gatk, dry_run=dry_run)
        if rc != 0:
            msg = f"GATK failed for {sample_name}: rc={rc} stderr={err.strip()}"
            logging.error(msg)
            return {"sample": sample_name, "ok": False, "message": msg}
    else:
        samtools_bin = which_or_path(samtools_exec or "samtools")
        if not samtools_bin:
            return {"sample": sample_name, "ok": False, "message": "samtools not found"}
        # samtools fastq -1 file1 -2 file2 -0 /dev/null -s /dev/null -n in.bam
        cmd_sam = [samtools_bin, "fastq", "-1", fastq_1, "-2", fastq_2, "-0", os.devnull, "-s", os.devnull, "-n", bam_path]
        rc, out, err = run_cmd(cmd_sam, dry_run=dry_run)
        if rc != 0:
            msg = f"samtools fastq failed for {sample_name}: rc={rc} stderr={err.strip()}"
            logging.error(msg)
            return {"sample": sample_name, "ok": False, "message": msg}

    # Check FASTQ files exist (unless dry-run)
    if not dry_run and (not os.path.exists(fastq_1) or not os.path.exists(fastq_2)):
        msg = f"FASTQ files not found for {sample_name} after conversion"
        logging.error(msg)
        return {"sample": sample_name, "ok": False, "message": msg}

    # 2) Run kallisto quant
    kallisto_bin = which_or_path(kallisto_exec or "kallisto")
    if not kallisto_bin:
        return {"sample": sample_name, "ok": False, "message": "kallisto not found"}

    cmd_kallisto = [
        kallisto_bin, "quant",
        "-i", kallisto_index,
        "-o", kallisto_sample_out,
        "--threads", str(kallisto_threads),
        fastq_1, fastq_2,
    ]

    rc, out, err = run_cmd(cmd_kallisto, dry_run=dry_run)
    if rc != 0:
        msg = f"kallisto failed for {sample_name}: rc={rc} stderr={err.strip()}"
        logging.error(msg)
        return {"sample": sample_name, "ok": False, "message": msg}

    # Optionally remove FASTQ files to save space
    if not keep_fastq and not dry_run:
        try:
            os.remove(fastq_1)
            os.remove(fastq_2)
        except Exception as e:
            logging.warning("Failed to delete FASTQ for %s: %s", sample_name, str(e))

    logging.info("Completed sample %s", sample_name)
    return {"sample": sample_name, "ok": True, "message": "Success"}


# --------------------------- CLI and orchestration -------------------------

def parse_args():
    p = argparse.ArgumentParser(description="BAM -> FASTQ -> Kallisto quantification pipeline")

    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam-list", help="Path to a text file listing BAM file paths (one per line).")
    group.add_argument("--bam-dir", help="Directory to scan recursively for .bam files.")

    p.add_argument("--fastq-out", required=True, help="Directory where FASTQ files will be written.")
    p.add_argument("--kallisto-out", required=True, help="Base directory for kallisto output (one subdir per sample).")
    p.add_argument("--kallisto-index", required=True, help="Path to kallisto transcriptome index file.")

    p.add_argument("--num-processes", type=int, default=4, help="Number of parallel worker processes.")
    p.add_argument("--kallisto-threads", type=int, default=4, help="Threads passed to kallisto per sample.")

    p.add_argument("--use-gatk", action="store_true", help="Use GATK SamToFastq (default: use samtools if --use-gatk not set).")
    p.add_argument("--gatk-exec", default="gatk", help="GATK executable name or full path.")
    p.add_argument("--samtools-exec", default="samtools", help="samtools executable name or full path.")
    p.add_argument("--kallisto-exec", default="kallisto", help="kallisto executable name or full path.")

    p.add_argument("--keep-fastq", action="store_true", help="Keep intermediate FASTQ files (do not delete them).")
    p.add_argument("--dry-run", action="store_true", help="Show commands without executing them.")
    p.add_argument("--log-file", default=None, help="Optional path to write log file.")

    return p.parse_args()


def discover_bams_from_dir(bam_dir: str) -> list[str]:
    bam_dir = Path(bam_dir)
    bam_files = [str(p) for p in bam_dir.rglob("*.bam")]
    return sorted(bam_files)


def main():
    args = parse_args()
    setup_logger(args.log_file)

    # Gather BAM files
    if args.bam_list:
        if not os.path.exists(args.bam_list):
            logging.error("bam-list file does not exist: %s", args.bam_list)
            sys.exit(1)
        with open(args.bam_list, "r") as fh:
            bam_files = [line.strip() for line in fh if line.strip()]
    else:
        bam_files = discover_bams_from_dir(args.bam_dir)

    if len(bam_files) == 0:
        logging.error("No BAM files found to process.")
        sys.exit(1)

    # Validate kallisto index
    if not os.path.exists(args.kallisto_index):
        logging.error("Kallisto index not found: %s", args.kallisto_index)
        sys.exit(1)

    # Make output directories
    os.makedirs(args.fastq_out, exist_ok=True)
    os.makedirs(args.kallisto_out, exist_ok=True)

    # Run pipeline in parallel
    results = []
    logging.info("Starting processing of %d BAM files with %d workers", len(bam_files), args.num_processes)

    with ProcessPoolExecutor(max_workers=args.num_processes) as exe:
        futures = {
            exe.submit(
                bam_to_fastq_and_kallisto,
                bam_path,
                args.fastq_out,
                args.kallisto_out,
                args.kallisto_index,
                args.use_gatk,
                args.gatk_exec,
                args.samtools_exec,
                args.kallisto_exec,
                args.kallisto_threads,
                args.keep_fastq,
                args.dry_run,
            ): bam_path for bam_path in bam_files
        }

        for fut in as_completed(futures):
            bam = futures[fut]
            try:
                res = fut.result()
                results.append(res)
                if res.get("ok"):
                    logging.info("%s: OK", res.get("sample"))
                else:
                    logging.error("%s: FAILED - %s", res.get("sample"), res.get("message"))
            except Exception as exc:
                logging.exception("Error processing %s: %s", bam, str(exc))

    # Summary
    ok_count = sum(1 for r in results if r.get("ok"))
    fail_count = len(results) - ok_count
    logging.info("Completed. Success: %d  Failures: %d", ok_count, fail_count)


if __name__ == "__main__":
    main()
