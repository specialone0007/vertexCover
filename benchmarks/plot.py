"""Render the benchmark CSVs and the quality table to PNG figures.

Usage:
    ./build/vertex_cover bench 50 --out benchmarks/results
    ./build/vertex_cover quality 500 > benchmarks/results/quality-500.csv
    python benchmarks/plot.py            # writes docs/figures/*.png

Only depends on matplotlib.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmarks" / "results"
FIGURES = ROOT / "docs" / "figures"

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
})
BLUE, RED, GREY = "#1f4e79", "#c0392b", "#7f8c8d"


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def plot_timing(path: Path, xlabel: str, title: str, out: Path) -> None:
    rows = read_csv(path)
    x = [r["size"] for r in rows]
    mean = [r["mean_ms"] for r in rows]
    lo = [r["ci95_low"] for r in rows]
    hi = [r["ci95_high"] for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.fill_between(x, lo, hi, color=BLUE, alpha=0.15, label="95 % CI")
    ax.plot(x, mean, "o-", color=BLUE, ms=4, label="mean")
    ax.set_xlabel(xlabel)
    ax.ticklabel_format(axis="x", style="plain")
    ax.set_ylabel("time per greedy run (ms)")
    ax.set_title(title, loc="left")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_quality(path: Path, out: Path) -> None:
    rows = read_csv(path)
    v = [r["V"] for r in rows]
    q = [r["mean_quality(opt/greedy)"] for r in rows]
    hits = [r["optimal_hits"] / r["trials"] * 100 for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(v, q, "o-", color=BLUE, ms=4, label="mean quality |OPT| / |greedy|")
    ax.axhline(1.0, color=GREY, lw=0.8, ls="--")
    ax.axhline(0.5, color=RED, lw=0.8, ls="--")
    ax.text(v[-1], 0.515, "2-approximation bound", ha="right", color=RED, fontsize=8)
    ax.set_ylim(0.45, 1.05)
    ax.set_xlabel("vertices V")
    ax.set_ylabel("quality")
    ax2 = ax.twinx()
    ax2.bar(v, hits, width=0.6, color=GREY, alpha=0.35, label="runs where greedy = OPT (%)")
    ax2.set_ylabel("greedy hits optimum (%)")
    ax2.set_ylim(0, 100)
    ax2.grid(False)
    ax2.spines["top"].set_visible(False)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper right")
    ax.set_title("Greedy vs exact on random graphs (500 per size)", loc="left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    timing = {
        "edges-fixed": ("vertices V  (E = 20,000)", "Running time, E fixed"),
        "vertices-fixed": ("edges E  (V = 20,000)", "Running time, V fixed"),
    }
    made = 0
    for prefix, (xlabel, title) in timing.items():
        files = sorted(RESULTS.glob(f"{prefix}-*-iter.csv"))
        if not files:
            print(f"no {prefix} CSV in {RESULTS}", file=sys.stderr)
            continue
        plot_timing(files[-1], xlabel, title, FIGURES / f"time-{prefix}.png")
        made += 1
    quality = sorted(RESULTS.glob("quality-*.csv"))
    if quality:
        plot_quality(quality[-1], FIGURES / "quality.png")
        made += 1
    print(f"wrote {made} figure(s) to {FIGURES}")
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())
