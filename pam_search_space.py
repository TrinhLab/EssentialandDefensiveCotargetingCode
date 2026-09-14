#!/usr/bin/env python3
"""
pam_search_space.py
====================
Counts non-target PAM search space for each gRNA in gRNA_Activity_Summary.xlsx
against a user-supplied S. aureus Newman genome FASTA.

Usage:
    python3 pam_search_space.py <genome.fasta> [options]

Arguments:
    genome.fasta        Path to S. aureus Newman genome in FASTA format
                        (single or multi-chromosome; all contigs are used)

Options:
    --xlsx PATH         Path to gRNA_Activity_Summary.xlsx
                        (default: same directory as this script)
    --seed-len INT      PAM-proximal seed length in bp (default: 5)
    --max-mm INT        Max mismatches allowed in seed match (default: 0)
    --out CSV           Output CSV filename (default: search_space_results.csv)
    --no-plot           Skip scatter plot

Background:
    SpCas9 recognises NGG PAM 3' of the protospacer on the non-target strand.
    Target interrogation initiates from the PAM-proximal (3') end of the spacer.
    The search space for a given gRNA is the count of genomic NGG sites where
    the immediately upstream (5') sequence matches the gRNA's PAM-proximal seed
    with <= max_mm mismatches. More matches = longer search time = potentially
    lower guide activity independent of gene function.
    (Moreb & Lynch, Nat. Commun. 2021; Sternberg et al., Nature 2014)
"""

import argparse
import csv
import sys
import os
import openpyxl

# ── Helpers ───────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def load_genome(fasta_path):
    """Load all contigs from FASTA into one concatenated uppercase string."""
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                continue
            seq_parts.append(line.upper())
    genome = ''.join(seq_parts)
    print(f"  Genome loaded: {len(genome):,} bp", flush=True)
    return genome

def load_grnas(xlsx_path):
    """
    Read Sheet1 from the xlsx.
    Expects columns: Gene, gRNA_ID, Sequence (20 nt spacer),
                     Sequence_and_PAM (23 nt), Log10R
    Returns list of dicts.
    """
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb['Sheet1']
    rows = list(ws.iter_rows(values_only=True))

    header_idx = next(i for i, r in enumerate(rows) if r[0] is not None)
    header = rows[header_idx]
    col = {h: i for i, h in enumerate(header) if h is not None}

    grnas = []
    for row in rows[header_idx + 1:]:
        if row[col['Gene']] is None:
            continue
        spacer = row[col['Sequence']]
        seq_pam = row[col['Sequence_and_PAM']]
        if not spacer or not seq_pam:
            continue
        grnas.append({
            'gene':    row[col['Gene']],
            'grna_id': row[col['gRNA_ID']],
            'spacer':  spacer.upper().strip(),
            'pam':     seq_pam[20:23].upper().strip(),
            'log10r':  row[col['Log10R']],
        })
    return grnas

def count_seed_matches(genome, seed_query, max_mm):
    """
    Count genomic NGG sites (both strands) where the seed_query matches
    the PAM-proximal protospacer sequence within max_mm mismatches.

    Forward strand pattern:  5'-[seed_query][NGG]-3'
    Reverse strand pattern:  5'-[CCN][RC(seed_query)]-3'  (CCN = NGG on rev strand)

    seed_query is the 3' end of the gRNA spacer (5'->3'), i.e. the last
    seed_len bases of the 20-nt spacer.
    """
    seed_len = len(seed_query)
    genome_len = len(genome)
    count = 0

    rc_seed = reverse_complement(seed_query)

    for i in range(genome_len - seed_len - 2):
        # Forward strand: [seed][xGG]
        pam_fwd = genome[i + seed_len : i + seed_len + 3]
        if pam_fwd[1] == 'G' and pam_fwd[2] == 'G':
            genomic_seed = genome[i : i + seed_len]
            if hamming(genomic_seed, seed_query) <= max_mm:
                count += 1

        # Reverse strand: [CCN][rc_seed]
        pam_rev = genome[i : i + 3]
        if pam_rev[0] == 'C' and pam_rev[1] == 'C':
            genomic_seed = genome[i + 3 : i + 3 + seed_len]
            if len(genomic_seed) == seed_len and hamming(genomic_seed, rc_seed) <= max_mm:
                count += 1

    return count

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Count PAM-proximal non-target search space for SpCas9 gRNAs.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('genome_pos', nargs='?', default=None,
                        help='Path to genome FASTA (positional)')
    parser.add_argument('-g', '--genome', default=None,
                        help='Path to genome FASTA')
    parser.add_argument('-x', '--xlsx', '--excel', default=None,
                        help='Path to gRNA_Activity_Summary.xlsx')
    parser.add_argument('-s', '--seed-len', '--seed', type=int, default=5,
                        help='PAM-proximal seed length in bp (default: 5)')
    parser.add_argument('-m', '--max-mm', '--mismatches', type=int, default=0,
                        help='Max mismatches in seed (default: 0)')
    parser.add_argument('-o', '--out', default='search_space_results.csv',
                        help='Output CSV (default: search_space_results.csv)')
    parser.add_argument('--no-plot', dest='plot', action='store_false',
                        default=True, help='Skip scatter plot')
    args = parser.parse_args()

    genome_path = args.genome or args.genome_pos
    if not genome_path:
        parser.error('Genome FASTA path is required (pass as first argument or via --genome / -g).')

    # Resolve xlsx
    if args.xlsx is None:
        candidates = [
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'gRNA_Activity_Summary.xlsx'),
            os.path.join(os.getcwd(), 'gRNA_Activity_Summary.xlsx'),
        ]
        for c in candidates:
            if os.path.exists(c):
                args.xlsx = c
                break
        if args.xlsx is None:
            sys.exit('ERROR: Cannot find gRNA_Activity_Summary.xlsx — pass --xlsx PATH')

    seed_len = args.seed_len
    max_mm   = args.max_mm

    print(f"Settings: seed_len={seed_len} bp, max_mismatches={max_mm}")
    print(f"gRNA file: {args.xlsx}")
    print(f"Genome:    {genome_path}\n")

    print("Loading genome...", flush=True)
    genome = load_genome(genome_path)

    print("Loading gRNAs...", flush=True)
    grnas = load_grnas(args.xlsx)
    print(f"  Loaded {len(grnas)} gRNAs\n")

    # Count seed matches per gRNA
    results = []
    for g in grnas:
        spacer = g['spacer']
        if len(spacer) < seed_len:
            print(f"  WARNING: {g['grna_id']} spacer too short, skipping")
            continue

        seed = spacer[-seed_len:]   # PAM-proximal = 3' end of spacer

        print(f"  {g['grna_id']:<14} seed: {seed}", end='', flush=True)
        n = count_seed_matches(genome, seed, max_mm)
        print(f"  → {n:5d} matches")

        results.append({
            'gene':        g['gene'],
            'grna_id':     g['grna_id'],
            'spacer':      spacer,
            'pam':         g['pam'],
            f'seed_{seed_len}bp': seed,
            'seed_matches': n,
            'log10r':       g['log10r'],
        })

    results.sort(key=lambda r: r['seed_matches'])

    # Write CSV
    fieldnames = ['gene', 'grna_id', 'spacer', 'pam',
                  f'seed_{seed_len}bp', 'seed_matches', 'log10r']
    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults written to: {args.out}")

    # Print summary table
    print(f"\n{'gRNA_ID':<14} {'Gene':<8} {'Seed':>{seed_len}} "
          f"{'Matches':>10}  {'Log10R':>8}")
    print('-' * (14 + 8 + seed_len + 10 + 8 + 6))
    for r in results:
        lr = f"{r['log10r']:.3f}" if r['log10r'] is not None else 'N/A'
        print(f"{r['grna_id']:<14} {r['gene']:<8} "
              f"{r[f'seed_{seed_len}bp']:>{seed_len}} "
              f"{r['seed_matches']:>10}  {lr:>8}")

    # Pearson r
    vals = [(r['seed_matches'], r['log10r'])
            for r in results if r['log10r'] is not None]
    r_pearson = None
    if len(vals) >= 3:
        import statistics
        xs = [v[0] for v in vals]
        ys = [v[1] for v in vals]
        xm, ym = statistics.mean(xs), statistics.mean(ys)
        sx = statistics.stdev(xs)
        sy = statistics.stdev(ys)
        if sx > 0 and sy > 0:
            r_pearson = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / \
                        ((len(xs) - 1) * sx * sy)
            print(f"\nPearson r (seed_matches vs Log10R): {r_pearson:+.3f}")
            if abs(r_pearson) < 0.3:
                interp = "Weak — search space does not strongly predict killing efficiency"
            elif r_pearson < 0:
                interp = "Negative — more search space correlates with LOWER killing (supports Moreb model)"
            else:
                interp = "Positive — more search space correlates with HIGHER killing (unexpected)"
            print(f"  → {interp}")

    # Scatter plot
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import numpy as np

            plot_vals = [(r['seed_matches'], r['log10r'], r['grna_id'])
                         for r in results if r['log10r'] is not None]
            xs = [v[0] for v in plot_vals]
            ys = [v[1] for v in plot_vals]
            labels = [v[2] for v in plot_vals]

            fig, ax = plt.subplots(figsize=(7.5, 5.5))
            ax.scatter(xs, ys, s=75, color='#2166ac',
                       edgecolors='white', linewidths=0.6, zorder=3)

            for x, y, lbl in zip(xs, ys, labels):
                ax.annotate(lbl, (x, y), textcoords='offset points',
                            xytext=(5, 3), fontsize=7.5, color='#333333')

            if r_pearson is not None and len(xs) >= 3:
                coeffs = np.polyfit(xs, ys, 1)
                x_line = np.linspace(min(xs), max(xs), 200)
                ax.plot(x_line, np.polyval(coeffs, x_line),
                        '--', color='#d6604d', linewidth=1.5,
                        label=f'Linear fit  r = {r_pearson:+.2f}')
                ax.legend(fontsize=9, framealpha=0.9)

            ax.set_xlabel(
                f'Non-target PAM search space\n'
                f'(NGG sites with {seed_len}-bp seed match, ≤{max_mm} mismatches)',
                fontsize=10)
            ax.set_ylabel('Log₁₀ reduction in viability (Fig. 3C)', fontsize=10)
            ax.set_title('PAM-proximal search space vs. gRNA killing efficiency\n'
                         'S. aureus Newman genome', fontsize=11, pad=8)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            plot_path = args.out.replace('.csv', '_scatter.png')
            plt.tight_layout()
            plt.savefig(plot_path, dpi=200, bbox_inches='tight')
            print(f"Scatter plot saved to: {plot_path}")

        except ImportError:
            print("matplotlib not available — skipping plot.")

if __name__ == '__main__':
    main()
