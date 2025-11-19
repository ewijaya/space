#!/usr/bin/env python3
"""
SPACE - Detection of Generic Spaced Motifs
Main execution script (Python port of run_all.pl)
"""

import sys
import os
import subprocess
from pathlib import Path

def main():
    # Parse arguments with defaults
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'test.fasta'
    species = sys.argv[2] if len(sys.argv) > 2 else 'SC'

    # Get base filename without extension
    base = Path(input_file).stem
    scname = f"{base}.sc_out"

    # Run multiparam.py to generate candidates and score
    print(f"Running SPACE analysis on {input_file} (species: {species})")
    with open(scname, 'w') as outfile:
        result = subprocess.run(
            ['python3', 'multiparam.py', input_file, species],
            stdout=outfile,
            stderr=subprocess.PIPE,
            text=True
        )

        if result.returncode != 0:
            print(f"Error in multiparam.py: {result.stderr}", file=sys.stderr)
            sys.exit(1)

    # Run advisor.py to find best motifs
    result = subprocess.run(
        ['python3', 'advisor.py', scname],
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        print(f"Error in advisor.py: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Clean up intermediate file
    try:
        os.unlink(scname)
    except FileNotFoundError:
        pass

if __name__ == '__main__':
    main()
