# Self-contained demo: plant a motif, scan it, test enrichment, make plots.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}


def random_seq(length, rng):
    return "".join(rng.choice(list(BASES), size=length))


def revcomp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def encode(seq):
    return np.array([IDX.get(b, -1) for b in seq], dtype=int)


def sample_motif(pfm, rng):
    # pfm columns are base frequencies; draw one instance.
    out = []
    for j in range(pfm.shape[1]):
        out.append(rng.choice(list(BASES), p=pfm[:, j]))
    return "".join(out)


def to_log_odds(freq, background, pseudocount=0.01):
    freq = freq + pseudocount
    freq = freq / freq.sum(axis=0, keepdims=True)
    bg = np.array(background).reshape(4, 1)
    return np.log2(freq / bg)


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


def permutation_pvalue(seqs, lom, n_perm, rng):
    observed = enrichment_stat(seqs, lom)
    null = np.empty(n_perm)
    for k in range(n_perm):
        shuffled = [shuffle_seq(s, rng) for s in seqs]
        null[k] = enrichment_stat(shuffled, lom)
    pval = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return observed, null, pval


def information_content(pfm):
    # Bits per position for a sequence logo.
    pfm = pfm + 1e-9
    pfm = pfm / pfm.sum(axis=0, keepdims=True)
    entropy = -np.sum(pfm * np.log2(pfm), axis=0)
    total_ic = np.log2(4) - entropy
    return pfm * total_ic


def plot_logo(ax, pfm):
    heights = information_content(pfm)
    colours = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}
    for j in range(heights.shape[1]):
        order = np.argsort(heights[:, j])
        bottom = 0.0
        for i in order:
            h = heights[i, j]
            if h <= 0:
                continue
            ax.text(j + 0.5, bottom + h / 2.0, BASES[i],
                    ha="center", va="center",
                    fontsize=6 + 40 * h / 2.0,
                    color=colours[BASES[i]], weight="bold")
            bottom += h
    ax.set_xlim(0, heights.shape[1])
    ax.set_ylim(0, 2)
    ax.set_ylabel("bits")
    ax.set_xlabel("position")
    ax.set_title("Sequence logo of planted motif")
    ax.set_xticks(np.arange(heights.shape[1]) + 0.5)
    ax.set_xticklabels(range(1, heights.shape[1] + 1))


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(42)

    # Define a motif PFM (A, C, G, T x width). A GC-rich-ish TF motif.
    pfm = np.array([
        [0.05, 0.80, 0.10, 0.05, 0.60, 0.10],
        [0.80, 0.05, 0.10, 0.05, 0.10, 0.70],
        [0.05, 0.10, 0.70, 0.05, 0.20, 0.10],
        [0.10, 0.05, 0.10, 0.85, 0.10, 0.10],
    ])
    pfm = pfm / pfm.sum(axis=0, keepdims=True)

    # Simulate promoter sequences; plant the motif in a fraction of them.
    n_seq = 40
    seq_len = 120
    plant_frac = 0.7
    seqs = []
    planted_flags = []
    for _ in range(n_seq):
        s = list(random_seq(seq_len, rng))
        if rng.random() < plant_frac:
            motif = sample_motif(pfm, rng)
            pos = rng.integers(0, seq_len - len(motif))
            s[pos:pos + len(motif)] = list(motif)
            planted_flags.append(True)
        else:
            planted_flags.append(False)
        seqs.append("".join(s))

    background = [0.25, 0.25, 0.25, 0.25]
    lom = to_log_odds(pfm, background)

    best = [best_score(s, lom) for s in seqs]
    observed, null, pval = permutation_pvalue(seqs, lom, 500, rng)

    # Save summary table.
    df = pd.DataFrame({
        "sequence_id": ["seq%d" % i for i in range(n_seq)],
        "planted": planted_flags,
        "best_score": np.round(best, 4),
    })
    df.to_csv("results/summary.csv", index=False)

    stats = pd.DataFrame({
        "metric": ["observed_stat", "null_mean", "permutation_pvalue"],
        "value": [round(observed, 4), round(float(np.mean(null)), 4), round(pval, 5)],
    })
    stats.to_csv("results/enrichment.csv", index=False)

    # Plot logo and null distribution.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    plot_logo(axes[0], pfm)

    axes[1].hist(null, bins=25, color="#888888", alpha=0.8, label="shuffled null")
    axes[1].axvline(observed, color="#d62728", lw=2,
                    label="observed (p=%.4f)" % pval)
    axes[1].set_xlabel("mean best PWM score")
    axes[1].set_ylabel("count")
    axes[1].set_title("Permutation enrichment test")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print("observed:", round(observed, 4), "p-value:", round(pval, 5))
    print("wrote figures/demo.png, results/summary.csv, results/enrichment.csv")


if __name__ == "__main__":
    main()
