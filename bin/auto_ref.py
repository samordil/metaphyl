#!/usr/bin/env python3
"""
Enhanced automatic reference genome selection with:
- Better reference selection metrics
- Improved error handling and logging
- Resource management
- More informative outputs
"""

import argparse
import datetime
import json
import logging
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

from Bio import SeqIO

# Constants
MINIMAP2_PRESET = "map-ont"  # Default for Oxford Nanopore reads
DEFAULT_MIN_COVERAGE = 80.0  # Minimum coverage percentage to consider
MAX_TMP_RETRIES = 3  # Max retries for temp file operations

class PipelineError(Exception):
    """Custom exception for pipeline failures with error codes"""
    def __init__(self, message: str, exit_code: int = 1):
        super().__init__(message)
        self.exit_code = exit_code

def configure_logging(log_file: str = "auto_ref.log") -> None:
    """Configure logging with both console and file output"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file)
        ]
    )
    logging.captureWarnings(True)  # Capture warnings as logs

def run_command(cmd: str, retries: int = 1) -> subprocess.CompletedProcess:
    """Execute shell command with robust error handling and retries"""
    for attempt in range(retries + 1):
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                executable='/bin/bash'
            )
            logging.debug(f"Command succeeded: {cmd}")
            return result
        except subprocess.CalledProcessError as e:
            if attempt == retries:
                logging.error(f"Command failed after {retries + 1} attempts: {cmd}")
                logging.error(f"Error output:\n{e.stderr}")
                raise PipelineError(f"Process failed: {e.stderr}")
            logging.warning(f"Command failed (attempt {attempt + 1}), retrying: {cmd}")
            continue

def calculate_reference_score(result: Dict[str, Any]) -> float:
    """Calculate a weighted score for reference selection"""
    weights = {
        'coverage': 0.5,    # Most important - percentage covered
        'meandepth': 0.3,   # Important - sequencing depth
        'numreads': 0.2     # Less important - raw read count
    }
    return (result['coverage'] * weights['coverage'] +
            result['meandepth'] * weights['meandepth'] +
            result['numreads'] * weights['numreads'] / 1000)  # Scale down large numbers

def process_reference(args: Tuple[str, str, Path, str]) -> Dict[str, Any]:
    """Process a single reference genome with improved error handling"""
    ref_id, ref_seq, reads, preset = args
    temp_files = []
    
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as ref_file:
            temp_files.append(ref_file.name)
            ref_file.write(f">{ref_id}\n{ref_seq}\n")
            ref_file.flush()

            # 1. Create index
            mmi_file = f"{ref_file.name}.mmi"
            run_command(f"minimap2 -d {mmi_file} {ref_file.name}")
            temp_files.append(mmi_file)
            
            # 2. Map reads and create sorted BAM
            sam_file = f"{ref_file.name}.sam"
            bam_file = f"{ref_file.name}.bam"
            sorted_bam = f"{ref_file.name}.sorted.bam"
            
            run_command(f"minimap2 -ax {preset} {mmi_file} {reads} > {sam_file}")
            temp_files.append(sam_file)
            
            run_command(f"samtools view -b {sam_file} > {bam_file}")
            temp_files.append(bam_file)
            
            run_command(f"samtools sort {bam_file} -o {sorted_bam}")
            temp_files.append(sorted_bam)
            
            # 3. Calculate coverage
            result = run_command(f"samtools coverage {sorted_bam}")
            
            # Parse coverage output
            lines = result.stdout.strip().split('\n')
            if len(lines) < 2:
                raise PipelineError("Incomplete coverage data from samtools")
            
            headers = [h.lower() for h in lines[0].split('\t')]  # Normalize headers
            data = lines[1].split('\t')
            
            if len(headers) != len(data):
                raise PipelineError("Mismatched coverage data columns")
            
            cov_data = dict(zip(headers, data))
            
            return {
                'ref_id': ref_id,
                'numreads': int(float(cov_data.get('numreads', 0))),
                'covbases': int(float(cov_data.get('covbases', 0))),
                'coverage': float(cov_data.get('coverage', 0)),
                'meandepth': float(cov_data.get('meandepth', 0)),
                'status': 'success',
                'score': calculate_reference_score({
                    'coverage': float(cov_data.get('coverage', 0)),
                    'meandepth': float(cov_data.get('meandepth', 0)),
                    'numreads': int(float(cov_data.get('numreads', 0)))
                })
            }
            
    except Exception as e:
        logging.warning(f"Failed processing {ref_id}: {str(e)}", exc_info=True)
        return {
            'ref_id': ref_id,
            'status': 'failed',
            'error': str(e)
        }
    finally:
        # Clean up temp files
        for f in temp_files:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception as e:
                logging.warning(f"Could not remove temp file {f}: {str(e)}")

def auto_map(reads: Path, msa: Path, output: Path, 
            threads: int, min_coverage: float = DEFAULT_MIN_COVERAGE) -> None:
    """Core mapping logic with improved reference selection"""
    try:
        # 1. Parse MSA
        logging.info(f"Parsing MSA: {msa}")
        refs = [(record.id, str(record.seq)) for record in SeqIO.parse(msa, "fasta")]
        if not refs:
            raise PipelineError("No references found in MSA", exit_code=10)

        # 2. Process references with progress tracking
        logging.info(f"Processing {len(refs)} references with {threads} threads")
        successful_results = []
        
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(
                    process_reference, 
                    (ref_id, seq, reads, MINIMAP2_PRESET)
                ): ref_id for ref_id, seq in refs
            }
            
            for future in as_completed(futures):
                ref_id = futures[future]
                try:
                    result = future.result()
                    if result['status'] == 'success':
                        successful_results.append(result)
                        logging.info(f"Completed {ref_id}: "
                                   f"{result['coverage']:.1f}% coverage, "
                                   f"{result['meandepth']:.1f}x depth")
                    else:
                        logging.warning(f"Failed {ref_id}: {result.get('error', 'Unknown error')}")
                except Exception as e:
                    logging.error(f"Unexpected error processing {ref_id}: {str(e)}")

        # 3. Select best reference
        if not successful_results:
            raise PipelineError("All reference mappings failed", exit_code=20)
        
        # Filter by minimum coverage
        filtered_results = [r for r in successful_results if r['coverage'] >= min_coverage]
        if not filtered_results:
            raise PipelineError(f"No references met minimum coverage of {min_coverage}%", exit_code=30)
        
        # Select by weighted score
        best_ref = max(filtered_results, key=lambda x: x['score'])
        
        logging.info(
            f"Selected reference: {best_ref['ref_id']}\n"
            f"  Coverage: {best_ref['coverage']:.2f}%\n"
            f"  Mean Depth: {best_ref['meandepth']:.1f}x\n"
            f"  Mapped Reads: {best_ref['numreads']}\n"
            f"  Selection Score: {best_ref['score']:.2f}"
        )

        # 4. Atomic write outputs
        temp_out = output.with_suffix('.tmp')
        try:
            with open(temp_out, 'w') as f:
                seq = next(r[1] for r in refs if r[0] == best_ref['ref_id'])
                f.write(f">{best_ref['ref_id']}\n{seq.replace('-', '')}\n")
            
            # Final atomic move
            temp_out.replace(output)
            
            # Write comprehensive metadata
            metadata = {
                'best_reference': best_ref['ref_id'],
                'selection_metrics': {
                    'coverage': best_ref['coverage'],
                    'meandepth': best_ref['meandepth'],
                    'numreads': best_ref['numreads'],
                    'score': best_ref['score']
                },
                'timestamp': datetime.datetime.now().isoformat(),
                'parameters': {
                    'reads_file': str(reads),
                    'msa_file': str(msa),
                    'min_coverage': min_coverage,
                    'threads': threads
                },
                'all_results': successful_results,
                'pipeline_version': '1.1.0'
            }
            
            with open(f"{output}.metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
                
        except Exception as e:
            if temp_out.exists():
                temp_out.unlink()
            raise PipelineError(f"Failed to write output: {str(e)}", exit_code=40)

    except Exception as e:
        logging.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        if isinstance(e, PipelineError):
            raise e
        raise PipelineError(str(e), exit_code=50)

def main():
    configure_logging()
    
    parser = argparse.ArgumentParser(
        description="Enhanced automated reference selection for genomic pipelines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--reads", type=Path, required=True,
                      help="Input reads (FASTQ/FASTA, can be gzipped)")
    parser.add_argument("--msa", type=Path, required=True,
                      help="Multiple sequence alignment (FASTA)")
    parser.add_argument("--output", type=Path, required=True,
                      help="Output reference file (FASTA)")
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1),
                      help="Number of threads to use (capped at 8 by default)")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                      help="Minimum coverage percentage to consider")
    parser.add_argument("--preset", type=str, default=MINIMAP2_PRESET,
                      help="Minimap2 preset for mapping")
    
    args = parser.parse_args()

    try:
        # Validate inputs
        if not args.reads.exists():
            raise PipelineError(f"Reads file not found: {args.reads}", exit_code=60)
        if not args.msa.exists():
            raise PipelineError(f"MSA file not found: {args.msa}", exit_code=70)
        
        args.output.parent.mkdir(parents=True, exist_ok=True)
        
        if args.min_coverage < 0 or args.min_coverage > 100:
            raise PipelineError("Minimum coverage must be between 0 and 100", exit_code=80)

        auto_map(args.reads, args.msa, args.output, args.threads, args.min_coverage)
        logging.info("Pipeline completed successfully")
        sys.exit(0)
        
    except PipelineError as e:
        logging.error(f"Pipeline aborted: {str(e)}")
        sys.exit(e.exit_code)
    except Exception as e:
        logging.error(f"Unexpected error: {str(e)}", exc_info=True)
        sys.exit(99)

if __name__ == "__main__":
    main()