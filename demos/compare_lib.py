import os, glob, argparse
from astropy.table import Table
import matplotlib.pyplot as pl

pl.style.use("via")
pl.rcParams["font.size"] = 14

MARKERS = {"MILES": "o", "C3K_LR": "s", "C3K_HR": "^"}
COLORS = {"MIST": "k", "PDVA": "firebrick", "PRSC": "dodgerblue", "BPSS": "orange"}


def get_parser():
    parser = argparse.ArgumentParser(description="Run pyFSPS AFE tests")
    parser.add_argument(
        "--dir",
        type=str,
        default="./lib_tests/sps_home-aa77c3_pyfsps-42635d_pyfslib-aa77c3",
        help="Output directory for test results",
    )
    parser.add_argument(
        "--kind",
        type=str,
        default="constneb",
        choices=["const", "constneb", "tau1"],
        help="Kind of SFH to plot spectrum of.",
    )
    parser.add_argument(
        "--filters",
        type=str,
        nargs=2,
        default=["bessell_B", "bessell_V"],
        help="Two filters to plot color evolution for.",
    )

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    dirn = args.dir
    f1, f2 = args.filters

    tbls = glob.glob(f"{dirn}/color_evol_*.csv")
    tags = [
        os.path.basename(t).replace("color_evol_", "").replace(".csv", "") for t in tbls
    ]
    # tags = [f"{t}_sfh1" for t in tags]
    tags = sorted(tags)
    print(tags)
    kind = args.kind

    fig, ax = pl.subplots(figsize=(8, 6))
    sfig, sax = pl.subplots(figsize=(8, 6))
    for tag in tags:
        isoc, slib, afe = tag.split("+")
        marker = MARKERS.get(slib.upper(), "o")
        color = COLORS.get(isoc.upper(), "k")
        color = "cyan" if afe == "afe1" else color
        color = "magenta" if afe == "afe1_afe0.3" else color
        try:
            table = Table.read(f"{dirn}/color_evol_{tag}.csv", format="csv")
            ax.plot(
                table["log_age"],
                table[f1] - table[f2],
                f"-{marker}",
                color=color,
                mec=color,
                alpha=0.5,
                label=f"{tag}, [a/Fe]=0.0",
            )
            spec = Table.read(f"{dirn}/spectrum_{tag}_{kind}.csv", format="csv")
            sax.plot(spec["wave"], spec["spec"], "-", alpha=0.5, label=f"{tag}")
        except Exception as e:
            print(f"Error plotting {tag}: {e}")
            pass

    ax.set_xlabel("SSP Age (log yr)")
    ax.set_ylabel("Color (B-V)")
    ax.set_ylim(-0.6, 1.8)
    ax.set_xlim(4.8, 10.4)
    ax.grid(True)
    # split the legend into two columns
    ax.legend(loc="upper left", frameon=True, framealpha=0.5, fontsize="small", ncol=2)

    sax.set_title("Spectra for different tags: constant SFH")
    sax.set_xlabel("Wavelength (AA)")
    sax.set_ylabel("Flux")
    sax.set_xlim(3000, 10000)
    sax.set_yscale("log")
    sax.set_ylim(5e-5, 5e-3)
    sax.grid(True)
    sax.legend(loc="upper right", frameon=True, framealpha=0.5, fontsize="small")
    fig.savefig(f"{dirn}/color_evol.png", dpi=300)
    sfig.savefig(f"{dirn}/spectrum_{kind}.png", dpi=300)
