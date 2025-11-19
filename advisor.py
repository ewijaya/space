#!/usr/bin/env python3
"""
SPACE - Motif Advisor
Python port of advisor.pl
Analyzes and ranks motif candidates
"""

import sys
import re
import math
from collections import defaultdict
from typing import List, Dict, Tuple, Set

# Configuration
IN_LIM = 0.6
TOP = 10
SBL_LLIM = 2
SBL_ULIM = 5000

def change_motif_to_iupac(amb_str: str) -> str:
    """Convert motif string with brackets to IUPAC notation."""
    iupac = {
        '[A]': 'A', '[T]': 'T', '[C]': 'C', '[G]': 'G',
        '[AC]': 'M', '[AG]': 'R', '[AT]': 'W', '[CG]': 'S',
        '[CT]': 'Y', '[GT]': 'K', '[ACG]': 'V', '[ACT]': 'H',
        '[AGT]': 'D', '[CGT]': 'B', '[ACGT]': 'N'
    }

    # Find all patterns like [ACGT] or single bases
    parts = re.findall(r'\[.*?\]|.', amb_str)
    result = []

    for part in parts:
        if part.startswith('['):
            # Sort bases inside brackets
            bases = ''.join(sorted(part[1:-1]))
            key = f'[{bases}]'
            result.append(iupac.get(key, part))
        else:
            result.append(part)

    return ''.join(result)


def reverse_iupac(amb_str: str) -> str:
    """Convert IUPAC codes back to bracket notation."""
    iupac = {
        'M': '[AC]', 'R': '[AG]', 'W': '[AT]', 'S': '[CG]',
        'Y': '[CT]', 'K': '[GT]', 'V': '[ACG]', 'H': '[ACT]',
        'D': '[AGT]', 'B': '[CGT]', 'N': '[ACGT]'
    }

    result = []
    for char in amb_str:
        result.append(iupac.get(char, char))
    return ''.join(result)


def bitwise(string: str) -> str:
    """Convert DNA sequence to bitwise representation."""
    bits = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
    parts = re.findall(r'\[.*?\]|.', string)
    result = []

    for part in parts:
        char_val = 0
        bases = re.findall(r'[ACGT]', part)
        for base in bases:
            char_val |= bits.get(base, 0)
        result.append(chr(char_val))

    return ''.join(result)


def matches(source: str, target: str) -> int:
    """Count matching positions between two bitwise strings."""
    mmatch = sum(1 for s, t in zip(source, target) if ord(s) & ord(t) == 0)
    return len(source) - mmatch


def get_match_count(source_: str, target_: str) -> int:
    """Get match count between two sequences (N matches anything)."""
    source = source_.replace('N', '[ATCG]')
    target = target_.replace('N', '[ATCG]')

    sc = bitwise(source)
    tg = bitwise(target)
    return matches(sc, tg)


def get_match_count_iupac(source_: str, target_: str) -> int:
    """Get match count with IUPAC codes."""
    source = reverse_iupac(source_)
    target = reverse_iupac(target_)

    sc = bitwise(source)
    tg = bitwise(target)
    return matches(sc, tg)


def motif_similarity(s1: str, s2: str) -> float:
    """Calculate similarity between two motifs."""
    # Determine shorter and longer
    if len(s1) <= len(s2):
        source, dest = s1, s2
    else:
        source, dest = s2, s1

    best = -1
    maxlen = 0

    for bp in range(len(dest)):
        for bp2 in range(len(source)):
            length = len(source) - bp2
            if (len(dest) - bp) < (len(source) - bp2):
                length = len(dest) - bp

            sstr = dest[bp:bp + length]
            sstr2 = source[bp2:bp2 + length]
            mcount = get_match_count_iupac(sstr, sstr2)

            if mcount > best:
                diff = len(dest) - length
                diff2 = len(source) - length
                maxlen = diff + length + diff2
                best = mcount

    if maxlen == 0:
        return 0.0

    sim = best / maxlen
    return round(sim, 3)


def find_sibling(mt_we_len: str, motif_list: List[str]) -> List[str]:
    """Find similar motifs (siblings)."""
    parts = mt_we_len.split()
    if len(parts) < 1:
        return []

    mt = parts[0]
    siblings = []

    for ls in motif_list:
        ls_parts = ls.split(' - ')
        if len(ls_parts) < 3:
            continue

        ls_param = ls_parts[1]
        ls_mtf = ls_parts[2]

        mtf_parts = ls_mtf.split()
        if len(mtf_parts) < 1:
            continue

        lmt = mtf_parts[0]

        if ls_mtf == mt_we_len:
            continue

        if motif_similarity(mt, lmt) > IN_LIM:
            siblings.append(f"{ls_param} - {ls_mtf}")

    return siblings


def parse_to_aoa(fname: str) -> Dict[str, List[List[str]]]:
    """Parse input file to array of arrays structure."""
    bighash = {}
    current_param = None
    temp = []

    try:
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip()

                # Match TAG: ParamGroupX
                match = re.match(r'^TAG:\s+(ParamGroup\d+)', line)
                if match:
                    current_param = match.group(1)
                    bighash[current_param] = []
                    continue

                # Match MOTIF : [ATCG]+ score
                match = re.match(r'^MOTIF\s?:\s+([\[\]ATCG]+)\s+(.*)', line)
                if match:
                    motif = match.group(1)
                    rest = match.group(2)
                    temp = [motif, rest]
                    continue

                # Match instance lines (number,-number,sequence)
                if re.match(r'^\d+,-\d+,[ATCGN]+$', line):
                    temp.append(line)
                    continue

                # Match separator
                if line.startswith('='):
                    if current_param and temp:
                        bighash[current_param].append(temp[:])
                    temp = []

    except FileNotFoundError:
        print(f"Error: File {fname} not found", file=sys.stderr)
        sys.exit(1)

    return bighash


def conv_to_ranked_hash(bh: Dict[str, List[List[str]]]) -> Dict[str, List]:
    """Convert to ranked hash structure."""
    nh = {}
    done = set()

    for pg in sorted(bh.keys()):
        ar = bh[pg]
        if not ar:
            continue

        count = 0
        for elem in ar:
            count += 1

            if not elem:
                continue

            motif = elem[0]
            rest = elem[1:] if len(elem) > 1 else []

            dstr = f"{motif}-{rest[0]}" if rest else motif

            if dstr in done:
                continue

            key = f"{count} Rank - {pg} - {motif}"
            nh[key] = rest
            done.add(dstr)

    return nh


def main():
    if len(sys.argv) < 2:
        print("Usage: advisor.py <input_file>", file=sys.stderr)
        sys.exit(1)

    filename = sys.argv[1]

    # Parse input file
    bhash = parse_to_aoa(filename)
    hash_data = conv_to_ranked_hash(bhash)

    # Filter by rank and score
    filtered_hash = {}
    for k, v in hash_data.items():
        parts = k.split(' - ')
        if len(parts) < 1:
            continue

        rank_str = parts[0]
        rank_match = re.search(r'(\d+)\s+Rank', rank_str)

        if rank_match:
            rankno = int(rank_match.group(1))
            score = float(v[0]) if v and len(v) > 0 else 0

            if rankno <= TOP and score > 0:
                filtered_hash[k] = v

    # Sort by score (descending), then rank (ascending)
    def sort_key(item):
        k, v = item
        score = float(v[0]) if v and len(v) > 0 else 0
        rank_match = re.search(r'(\d+)\s+Rank', k)
        rank = int(rank_match.group(1)) if rank_match else 999
        return (-score, rank)

    motif_set = sorted(filtered_hash.items(), key=sort_key)
    motif_set_keys = [k for k, v in motif_set]

    # Find siblings
    sibling_set = {}
    for ms in motif_set_keys:
        parts = ms.split(' - ')
        if len(parts) >= 3:
            motif = parts[2]
            sbls = find_sibling(motif, motif_set_keys)
            sibling_set[ms] = sbls

    # Print results
    cnt = 0
    for mset in motif_set_keys:
        no_of_sbl = len(sibling_set.get(mset, []))

        if SBL_LLIM <= no_of_sbl < SBL_ULIM:
            cnt += 1

            parts = mset.split(' - ')
            if len(parts) < 3:
                continue

            mot_e_len = parts[2]
            mot_parts = mot_e_len.split()

            if cnt == 1:
                print("==== BEST Motif Seems To Be: ===")
                if mot_parts:
                    mot = mot_parts[0]
                    print(f"MOTIF: {change_motif_to_iupac(mot)}\n")

                    wscr = filtered_hash[mset][0] if filtered_hash[mset] else 0
                    print(f"Significance Score : {wscr}")

                    if len(filtered_hash[mset]) > 1:
                        inset = filtered_hash[mset][1:]
                        print('\n'.join(inset))
                    print("\n")
            else:
                print("==== Other Potentially Good Motif: ===")
                if mot_parts:
                    mot = mot_parts[0]
                    print(f"MOTIF: {change_motif_to_iupac(mot)}\n")

                    if filtered_hash[mset]:
                        print(f"Significance Score : {filtered_hash[mset][0]}")
                        if len(filtered_hash[mset]) > 1:
                            print('\n'.join(filtered_hash[mset][1:]))
                    print()


if __name__ == '__main__':
    main()
