# Command-line PWM motif enrichment scanner.
import argparse
import numpy as np

BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}


def read_fasta(path):
    seqs = []
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(buf))
                name = line[1:]
                buf = []
            else:
                buf.append(line.upper())
    if name is not None:
        seqs.append("".join(buf))
    return seqs


def read_pwm(path):
    # Expect 4 whitespace-separated rows (A, C, G, T) of counts or frequencies.
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])
    mat = np.array(rows, dtype=float)
    if mat.shape[0] != 4:
        raise ValueError("PWM must have 4 rows in A, C, G, T order")
    return mat


def to_log_odds(freq, background, pseudocount=0.01):
    freq = freq + pseudocount
    freq = freq / freq.sum(axis=0, keepdims=True)
    bg = np.array(background).reshape(4, 1)
    return np.log2(freq / bg)


def encode(seq):
    return np.array([IDX.get(b, -1) for b in seq], dtype=int)


def revcomp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def best_score(seq, lom):
    w = lom.shape[1]
    best = -np.inf
    for strand in (seq, revcomp(seq)):
        enc = encode(strand)
        for i in range(len(enc) - w + 1):
            window = enc[i:i + w]
            if np.any(window < 0):
                continue
            s = lom[window, np.arange(w)].sum()
            if s > best:
                best = s
    return best


def enrichment_stat(seqs, lom):
    return float(np.mean([best_score(s, lom) for s in seqs]))


def shuffle_seq(seq, rng):
    arr = list(seq)
    rng.shuffle(arr)
    return "".join(arr)


def permutation_pvalue(seqs, lom, n_perm, seed=0):
    rng = np.random.default_rng(seed)
    observed = enrichment_stat(seqs, lom)
    null = np.empty(n_perm)
    for k in range(n_perm):
        shuffled = [shuffle_seq(s, rng) for s in seqs]
        null[k] = enrichment_stat(shuffled, lom)
    pval = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return observed, null, pval


def main():
    ap = argparse.ArgumentParser(description="PWM motif enrichment scanner")
    ap.add_argument("--fasta", required=True, help="target sequences (FASTA)")
    ap.add_argument("--pwm", required=True, help="PWM file, 4 rows A C G T")
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    pwm = read_pwm(args.pwm)
    background = [0.25, 0.25, 0.25, 0.25]
    lom = to_log_odds(pwm, background)
    observed, null, pval = permutation_pvalue(seqs, lom, args.perm, args.seed)
    print("sequences:", len(seqs))
    print("observed enrichment statistic:", round(observed, 4))
    print("null mean:", round(float(np.mean(null)), 4))
    print("permutation p-value:", round(pval, 5))


if __name__ == "__main__":
    main()
