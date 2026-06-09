import argparse
import itertools
import pathlib

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate plots for hmatrix build benchmark."
    )

    parser.add_argument(
        "--inputfile",
        type=str,
        help="Path to input file",
    )

    parser.add_argument(
        "--outputfile",
        type=str,
        help="Path to output file",
    )

    parser.add_argument(
        "--show", action="store_true", help="Display the plot interactively"
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputfile = pathlib.Path(args.inputfile)
    outputfile = pathlib.Path(args.outputfile)
    show = args.show

    data = pd.read_csv(inputfile)
    data.columns = data.columns.str.replace(" ", "")

    # Parameters
    epsilons = data["epsilon"].unique()
    generators = data["generator_type"].unique()
    policies = data["policy_type"].unique()
    sizes = data["size"].unique()
    low_rank_generators = data["low_rank_generator_type"].unique()
    clustering_types = data["clustering_type"].unique()

    markers = itertools.cycle(["o", "s", "^", "v", "x", "d", "*"])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r"Size $N$")
    ax.set_ylabel("Time (s)")
    ax.set_title(r"Time for $\mathcal{H}$matrix assembly with sizes " + str(sizes))
    for epsilon in epsilons:
        for generator_type in generators:
            for policy in policies:
                for low_rank_generator in low_rank_generators:
                    for clustering_type in clustering_types:
                        local_data = data[
                            (data["generator_type"] == generator_type)
                            & (data["epsilon"] == epsilon)
                            & (data["id_rep"] == "mean")
                            & (data["policy_type"] == policy)
                            & (data["low_rank_generator_type"] == low_rank_generator)
                            & (data["clustering_type"] == clustering_type)
                        ]
                        ax.loglog(
                            local_data["size"],
                            local_data["time(s)"],
                            marker=next(markers),
                            label=r"$\varepsilon=$"
                            + f"{epsilon}, {generator_type}, {policy}, {low_rank_generator}, {clustering_type}",
                        )

    ref_data = data[
        (data["generator_type"] == generators[0])
        & (data["epsilon"] == epsilons[0])
        & (data["policy_type"] == policies[0])
        & (data["id_rep"] == "mean")
    ]
    ax.plot(
        ref_data["size"],
        np.log(ref_data["size"])
        * np.log(ref_data["size"])
        * ref_data["size"]
        * ref_data["time(s)"].iloc[0]
        / (
            np.log(ref_data["size"].iloc[0])
            * np.log(ref_data["size"].iloc[0])
            * ref_data["size"].iloc[0]
        ),
        label=r"$N \log^2 N$",
        color="black",
        linestyle=":",
    )
    ax.legend()

    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())

    ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

    ax = fig.add_subplot(211)
    ax.set_xlabel(r"Size $N$")
    ax.set_ylabel("Compression")
    ax.set_title(r"Compression for $\mathcal{H}$matrix with sizes " + str(sizes))
    for epsilon in epsilons:
        for generator_type in generators:
            for policy in policies:
                for low_rank_generator in low_rank_generators:
                    for clustering_type in clustering_types:
                        local_data = data[
                            (data["generator_type"] == generator_type)
                            & (data["epsilon"] == epsilon)
                            & (data["id_rep"] == "mean")
                            & (data["policy_type"] == policy)
                            & (data["low_rank_generator_type"] == low_rank_generator)
                            & (data["clustering_type"] == clustering_type)
                        ]
                        ax.plot(
                            local_data["size"],
                            local_data["space_saving"],
                            marker=next(markers),
                            label=r"$\varepsilon=$"
                            + f"{epsilon}, {generator_type}, {policy}, {low_rank_generator}, {clustering_type}",
                        )

    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())

    ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

    if show:
        plt.show()

    fig.savefig(outputfile)


if __name__ == "__main__":
    main()
