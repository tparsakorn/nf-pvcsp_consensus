#!/usr/bin/env python3
"""
check_orf.py — QC gate for PvCSP NGSpeciesID consensus sequences.

Validation only: does not modify sequence.

The current pipeline does not strip PCR-primer sequence from the final
consensus (--primer_file is not set), so "position 0" of these FASTA
records is not the true CDS start, and orientation (which strand carries
the coding sequence) is not guaranteed either. A naive fixed-frame,
forward-strand-only stop-codon scan is therefore unreliable here — it will
pick up stop codons that fall in non-coding primer-derived flanking
sequence and misreport a clean haplotype as broken.

Instead this anchors on a short peptide motif that is conserved between
VK210 and VK247 immediately upstream of the central repeat
('RENKLKQP'), tries both strands and all 3 frames, and evaluates stop
codons / repeat-unit structure only in the sequence downstream of that
anchor — which is exactly the CDS region that matters for this assay.

Usage:
    python3 check_orf.py consensus_1.fasta consensus_2.fasta ...
    python3 check_orf.py --glob "barcode*_pvcsp/medaka_cl_id_*/consensus.fasta"
"""
import sys
import argparse
import glob
import re

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

ANCHOR = "RENKLKQP"   # conserved, immediately upstream of the central repeat in both types
VK210_UNIT = re.compile(r'GDRA[DA]GQPA')
VK247_UNIT = re.compile(r'ANGA[GD][NDG]QPG')


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def translate(seq, frame):
    seq = seq.upper().replace('U', 'T')
    return ''.join(CODON_TABLE.get(seq[i:i+3], 'X') for i in range(frame, len(seq) - 2, 3))


def read_fastx(path):
    recs = []
    with open(path) as fh:
        first = fh.readline()
        fh.seek(0)
        if first.startswith('@'):
            lines = [l.rstrip('\n') for l in fh]
            for i in range(0, len(lines), 4):
                recs.append((lines[i][1:].split()[0], lines[i+1]))
        else:
            name, seq = None, []
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if name is not None:
                        recs.append((name, ''.join(seq)))
                    name, seq = line[1:].split()[0], []
                else:
                    seq.append(line)
            if name is not None:
                recs.append((name, ''.join(seq)))
    return recs


def find_anchor(prot):
    idx = prot.find(ANCHOR)
    if idx != -1:
        return idx
    for i in range(max(0, len(prot) - len(ANCHOR) + 1)):
        window = prot[i:i + len(ANCHOR)]
        if sum(1 for a, b in zip(window, ANCHOR) if a != b) <= 1:
            return i
    return -1


def analyze(seq):
    """Try both strands / all 3 frames, anchor on the conserved motif,
    return the best (fewest post-anchor stop codons) candidate."""
    best = None
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in (0, 1, 2):
            prot = translate(s, frame)
            idx = find_anchor(prot)
            if idx == -1:
                continue
            after = prot[idx:]
            first_stop = after.find('*')
            # A stop within the last few residues is the true/terminal stop
            # followed by a little untrimmed 3' primer readthrough (this
            # pipeline does not run --primer_file), not an internal
            # frameshift. Only count stops that leave a meaningful amount
            # of sequence after them as "internal".
            TERMINAL_SLACK = 6
            n_stops = sum(1 for i, aa in enumerate(after)
                          if aa == '*' and i < len(after) - TERMINAL_SLACK)
            score = (n_stops, -len(after))
            cand = (score, strand, frame, len(after), n_stops, after)
            if best is None or cand < best:
                best = cand
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--glob', dest='glob_pat', default=None)
    args = ap.parse_args()

    files = list(args.files)
    if args.glob_pat:
        files += sorted(glob.glob(args.glob_pat, recursive=True))
    if not files:
        ap.print_help()
        sys.exit(1)

    header = f"{'file':<50}{'id':<24}{'strand':>7}{'frame':>6}{'post_aa':>8}{'stops':>6}{'rep_units':>10}  verdict"
    print(header)
    print('-' * len(header))

    any_fail = False
    for path in files:
        for name, seq in read_fastx(path):
            r = analyze(seq)
            if r is None:
                print(f"{path:<50}{name:<24}  NO ANCHOR FOUND (len={len(seq)}) -- primer/adapter "
                      f"may still be masking the anchor, or this is not a PvCSP consensus")
                any_fail = True
                continue
            score, strand, frame, alen, n_stops, after = r
            vk210 = len(VK210_UNIT.findall(after))
            vk247 = len(VK247_UNIT.findall(after))
            rep_units = f"VK210x{vk210}" if vk210 > vk247 else (f"VK247x{vk247}" if vk247 else "0")
            ok = n_stops == 0
            verdict = "OK" if ok else "FAIL"
            if not ok:
                any_fail = True
            print(f"{path:<50}{name:<24}{strand:>7}{frame:>6}{alen:>8}{n_stops:>6}{rep_units:>10}  {verdict}")

    sys.exit(1 if any_fail else 0)


if __name__ == '__main__':
    main()
