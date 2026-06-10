#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def read_output(path):
    i_vals = []
    rho = []
    vel = []
    pres = []
    rho_exact = []
    vel_exact = []
    pres_exact = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("i"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            i_vals.append(int(parts[0]))
            rho.append(float(parts[1]))
            vel.append(float(parts[2]))
            pres.append(float(parts[3]))
            if len(parts) >= 7:
                rho_exact.append(float(parts[4]))
                vel_exact.append(float(parts[5]))
                pres_exact.append(float(parts[6]))

    if not i_vals:
        raise ValueError(f"no data found in {path}")

    n_cell = max(i_vals) + 1
    x = [(idx + 0.5) / n_cell for idx in i_vals]

    exact = None
    if rho_exact:
        exact = (rho_exact, vel_exact, pres_exact)

    return x, rho, vel, pres, exact


def plot_output(path, output=None, show=True):
    x, rho, vel, pres, exact = read_output(path)

    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    axes[0].plot(x, rho, "o-", markersize=4, label="numerical")
    axes[0].set_ylabel("density")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(x, vel, "o-", markersize=4, color="tab:orange", label="numerical")
    axes[1].set_ylabel("velocity")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(x, pres, "o-", markersize=4, color="tab:green", label="numerical")
    axes[2].set_ylabel("pressure")
    axes[2].set_xlabel("x")
    axes[2].grid(True, alpha=0.3)

    if exact is not None:
        rho_exact, vel_exact, pres_exact = exact
        axes[0].plot(x, rho_exact, "-", linewidth=1.5, color="tab:red", label="exact")
        axes[1].plot(x, vel_exact, "-", linewidth=1.5, color="tab:red", label="exact")
        axes[2].plot(x, pres_exact, "-", linewidth=1.5, color="tab:red", label="exact")
        for ax in axes:
            ax.legend(loc="best")

    fig.suptitle(f"Sod shocktube: {path.name}")
    fig.tight_layout()

    if output is not None:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"wrote {output}")

    if show:
        plt.show()

    return fig


def main():
    parser = argparse.ArgumentParser(
        prog="plot_output.py",
        description="Plot density, velocity, and pressure from shocktube output.txt.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  plot_output.py
  plot_output.py output.txt
  plot_output.py -o sod.png --no-show
""",
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="output.txt",
        help="input data file (default: output.txt)",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="save figure to FILE",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="do not open an interactive plot window",
    )
    args = parser.parse_args()

    path = Path(args.input)
    plot_output(path, output=args.output, show=not args.no_show)


if __name__ == "__main__":
    main()
