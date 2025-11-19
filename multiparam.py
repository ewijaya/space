#!/usr/bin/env python3
"""
SPACE - Multiple Parameter Configuration Runner
Python port of multiparam.pl
"""

import sys
import os
import subprocess
from pathlib import Path
import math

def main():
    # Parse arguments with defaults
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'test.fasta'
    species = sys.argv[2] if len(sys.argv) > 2 else 'MM'

    submot_len = 5
    base = Path(input_file).stem

    # Parameter combinations to try
    cov_list = [0.5, 0.8]
    quorum_list = [0.5, 0.9]
    wlist = [10, 15]

    param_count = 0

    for w in wlist:
        for q in quorum_list:
            for cov in cov_list:
                param_count += 1

                ncov = math.floor(cov * w)

                # Create config file
                config_file = f"{base}-{w}-{submot_len}-{ncov}-{q}.config"
                print(f"{config_file} PARAM {param_count}")

                # Determine paths for config file
                input_file_for_config = input_file if not input_file.startswith('files/') else input_file.replace('files/', '', 1)

                with open(config_file, 'w') as f:
                    f.write(f"{submot_len}\n")
                    f.write(f"{ncov}\n")
                    f.write(f"{q}\n")
                    f.write(f"{w}\n")
                    f.write(f"FreqFiles/{species}.6.freq\n")
                    f.write(f"FreqFiles/{species}.8.freq\n")
                    f.write(f"files/{input_file_for_config}\n")
                    f.write(f"{base}.dat\n")

                # File paths
                base_out = f"{base}.out"
                base_dat = f"{base}.dat"
                # Add files/ prefix only if not already present
                if input_file.startswith('files/'):
                    dir_input_file = input_file
                else:
                    dir_input_file = f"files/{input_file}"

                # Run generate_v.exe
                with open(dir_input_file, 'r') as infile, open(base_out, 'w') as outfile:
                    subprocess.run(
                        ['./generate_v.exe', config_file],
                        stdin=infile,
                        stdout=outfile,
                        check=True
                    )

                # Run mining_v.exe
                with open(base_out, 'r') as infile, open(base_dat, 'w') as outfile:
                    subprocess.run(
                        ['./mining_v.exe'],
                        stdin=infile,
                        stdout=outfile,
                        check=True
                    )

                # Run scoring_v.exe
                print(f"TAG: ParamGroup{param_count}")
                subprocess.run(
                    ['./scoring_v.exe', config_file],
                    check=True
                )

                print("\n")

                # Clean up temporary files
                for temp_file in [config_file, base_out, base_dat]:
                    try:
                        os.unlink(temp_file)
                    except FileNotFoundError:
                        pass

if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
