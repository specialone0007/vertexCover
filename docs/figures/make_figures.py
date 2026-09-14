"""Static illustrations for the report and README (no networkx needed).

    python docs/figures/make_figures.py
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parent
BLUE, RED, GREEN, GREY, INK = "#1f4e79", "#c0392b", "#27ae60", "#bdc3c7", "#2c3e50"
plt.rcParams.update({"figure.dpi": 170, "font.size": 9})


def node(ax, xy, label, fill="white", edge=INK, r=0.22, fs=8, text=INK):
    ax.add_patch(plt.Circle(xy, r, facecolor=fill, edgecolor=edge, lw=1.3, zorder=3))
    ax.text(*xy, label, ha="center", va="center", fontsize=fs, color=text, zorder=4)


def edge(ax, a, b, color=INK, lw=1.1, ls="-"):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls, zorder=2, solid_capstyle="round")


def clean(ax):
    ax.set_aspect("equal")
    ax.axis("off")


# --------------------------------------------------------------------------- greedy vs optimum
def star_and_petersen() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.7))

    # Star K1,5 — optimum is the centre; greedy takes an edge
    ax = axes[0]
    c = (0, 0)
    leaves = [(math.cos(a) * 1.1, math.sin(a) * 1.1) for a in
              [math.pi / 2 + i * 2 * math.pi / 5 for i in range(5)]]
    for p in leaves:
        edge(ax, c, p)
    for i, p in enumerate(leaves):
        node(ax, p, str(i + 1), fill=RED if i == 0 else "white", edge=RED if i == 0 else INK,
             text="white" if i == 0 else INK)
    node(ax, c, "0", fill=RED, edge=RED, text="white")
    ax.set_title("Star $K_{1,5}$\ngreedy = 2, optimum = 1 (ratio 2)", fontsize=8.5)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); clean(ax)

    # Petersen graph — greedy 10 vs optimum 6
    def petersen(ax, cover, title):
        outer = [(math.cos(a) * 1.25, math.sin(a) * 1.25) for a in
                 [math.pi / 2 + i * 2 * math.pi / 5 for i in range(5)]]
        inner = [(math.cos(a) * 0.6, math.sin(a) * 0.6) for a in
                 [math.pi / 2 + i * 2 * math.pi / 5 for i in range(5)]]
        pos = outer + inner
        E = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
             (5, 7), (7, 9), (9, 6), (6, 8), (8, 5)]
        for u, v in E:
            edge(ax, pos[u], pos[v])
        for i, p in enumerate(pos):
            inc = i in cover
            node(ax, p, str(i), fill=RED if inc else "white", edge=RED if inc else INK,
                 text="white" if inc else INK, r=0.19, fs=7)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); clean(ax)

    petersen(axes[1], {0, 1, 2, 3, 4, 9, 5, 7, 6, 8}, "Petersen graph\ngreedy cover: 10 vertices")
    petersen(axes[2], {1, 3, 4, 5, 6, 7}, "Petersen graph\nminimum cover: 6 vertices")
    fig.tight_layout()
    fig.savefig(OUT / "greedy-vs-optimum.png")
    plt.close(fig)


# --------------------------------------------------------------------------- 3-SAT reduction gadget
def reduction_gadget() -> None:
    """(x1 ∨ x1 ∨ x2) ∧ (¬x1 ∨ ¬x2 ∨ ¬x2) ∧ (¬x1 ∨ x2 ∨ x2), k = L + 2M = 2 + 6 = 8.
    Satisfying assignment x1 = F, x2 = T.  Red = vertices in the cover."""
    fig, ax = plt.subplots(figsize=(8.2, 3.6))

    # literal pairs on top
    lit = {"x1": (1.0, 2.6), "!x1": (2.2, 2.6), "x2": (5.4, 2.6), "!x2": (6.6, 2.6)}
    truth = {"x1": False, "!x1": True, "x2": True, "!x2": False}
    edge(ax, lit["x1"], lit["!x1"]); edge(ax, lit["x2"], lit["!x2"])

    # gadgets (triangles) at bottom
    clauses = [("x1", "x1", "x2"), ("!x1", "!x2", "!x2"), ("!x1", "x2", "x2")]
    centres = [1.4, 3.8, 6.2]
    tri_pts = []
    for cx, cl in zip(centres, clauses):
        pts = [(cx - 0.55, 0.35), (cx + 0.55, 0.35), (cx, 1.25)]
        tri_pts.append((pts, cl))
        for i in range(3):
            edge(ax, pts[i], pts[(i + 1) % 3])
        for p, l in zip(pts, cl):
            edge(ax, p, lit[l], color=GREY, lw=0.9)

    # cover: true literals + the two triangle vertices NOT connected to a true literal
    # (if several qualify, pick any two).
    for name, p in lit.items():
        t = truth[name]
        node(ax, p, name.replace("!", "¬").replace("1", "₁").replace("2", "₂"), fill=RED if t else "white",
             edge=RED if t else INK, text="white" if t else INK, r=0.26)
    for pts, cl in tri_pts:
        satisfied = [i for i, l in enumerate(cl) if truth[l]]
        keep_out = satisfied[0]           # one vertex whose outside edge is already covered
        for i, (p, l) in enumerate(zip(pts, cl)):
            inc = i != keep_out
            node(ax, p, l.replace("!", "¬").replace("1", "₁").replace("2", "₂"), fill=RED if inc else "white",
                 edge=RED if inc else INK, text="white" if inc else INK, r=0.24, fs=7)

    ax.text(0.15, 3.1, "literal gadgets  (one of each pair must be in the cover)", fontsize=8, color=INK)
    ax.text(0.15, -0.35, "clause gadgets  (triangles: at least two of three vertices in any cover)",
            fontsize=8, color=INK)
    ax.text(8.0, 3.5, "φ = (x₁∨x₁∨x₂) ∧ (¬x₁∨¬x₂∨¬x₂) ∧ (¬x₁∨x₂∨x₂)\n"
                      "k = L + 2M = 2 + 2·3 = 8\n"
                      "x₁ = F, x₂ = T  →  cover of size 8 (red)",
            fontsize=8, color=INK, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f4f6f7", ec=GREY))
    ax.set_xlim(-0.1, 8.1); ax.set_ylim(-0.6, 3.6); clean(ax)
    fig.tight_layout()
    fig.savefig(OUT / "reduction-3sat.png")
    plt.close(fig)


if __name__ == "__main__":
    star_and_petersen()
    reduction_gadget()
    print("wrote greedy-vs-optimum.png, reduction-3sat.png")
