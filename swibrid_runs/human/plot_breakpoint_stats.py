

def setup_argparse(parser):
    parser.add_argument(
        "-g",
        "--gaps",
        dest="gaps",
        help="""required: file with gap positions (output of get_gaps.py)""",
    )
    parser.add_argument(
        "-c",
        "--clustering",
        dest="clustering",
        help="""file with clustering results""",
    )
    parser.add_argument(
        "-a",
        "--clustering_analysis",
        dest="clustering_analysis",
        help="""file with clustering analysis""",
    )
    parser.add_argument("-o", "--output", help="""output plot""")
    parser.add_argument(
        "-b",
        "--binsize",
        dest="binsize",
        default=50,
        type=int,
        help="""binsize [50]""",
    )
    parser.add_argument(
        "--scale_factor",
        dest="scale_factor",
        default=10,
        type=int,
        help="""factor to increase binsize for 2D histogram [10]""",
    )
    parser.add_argument(
        "--max_gap",
        dest="max_gap",
        default=75,
        type=int,
        help="""max gap size to ignore [75]""",
    )
    parser.add_argument("--sample", dest="sample", help="""sample name (for figure title)""")
    parser.add_argument(
        "--switch_coords",
        dest="switch_coords",
        default="chr14:106050000-106337000:-",
        help="""coordinates of switch region [chr14:106050000-106337000:-]""",
    )
    parser.add_argument(
        "--switch_annotation",
        dest="switch_annotation",
        help="""bed file with switch annotation""",
    )
    parser.add_argument(
        "--use_clones",
        dest="use_clones",
        help="""comma-separated list of clones to use or 'all' (default: filtered clusters, excluding singletons)""",
    )
    parser.add_argument(
        "--weights",
        dest="weights",
        default="cluster",
        help="""use different weights ("cluster" | "reads" | "adjusted") [cluster]""",
    )


if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import scipy.stats
    import pysam
    from logzero import logger
    from swibrid.scripts.utils import (
        parse_switch_coords,
        read_switch_anno,
        shift_coord,
        get_switch_coverage,
        parse_range,
        get_switch_iis,
    )
    from argparse import ArgumentParser

    parser=ArgumentParser()
    setup_argparse(parser)
    args=parser.parse_args()

    (
        switch_chrom,
        switch_start,
        switch_end,
        switch_orientation,
    ) = parse_switch_coords(args.switch_coords)
    switch_anno = read_switch_anno(args.switch_annotation)
    cov_int, Ltot, eff_start, eff_end, anno_recs = get_switch_coverage(
        switch_anno, switch_chrom, switch_start, switch_end
    )

    binsize = args.binsize
    switch_iis = get_switch_iis(anno_recs, cov_int, eff_start, binsize)

    logger.info("reading gaps from " + args.gaps)
    gaps = np.load(args.gaps)
    nreads = gaps["read_idx"].max() + 1

    if args.clustering:
        logger.info("reading clustering from " + args.clustering)
        clustering = pd.read_csv(args.clustering, header=0, index_col=0)
        clustering = clustering[~clustering["cluster"].isna()]
        nreads = clustering.shape[0]
        logger.info("reading clustering analysis from " + args.clustering_analysis)
        analysis = pd.read_csv(args.clustering_analysis, header=0, index_col=0)
        if args.use_clones:
            if args.use_clones == "all":
                clones = clustering["cluster"].astype(int).unique()
            else:
                clones = list(map(int, args.use_clones.split(",")))
        else:
            clusters = clustering["filtered_cluster"].dropna()
            clones = clusters[clusters >= 0].astype(int).unique()
        logger.info("using {0} clones".format(len(clones)))
        nc = clustering["cluster"].dropna().astype(int).value_counts()
        singletons = nc.index[nc == 1]
        if args.weights == "cluster":
            logger.info("using uniform weights per cluster")
            w = 1.0 / nc.loc[clustering["cluster"].values]
        elif args.weights == "reads":
            logger.info("using uniform weights per read")
            w = pd.Series(1, index=np.arange(nreads))
        elif args.weights == "adjusted":
            logger.info("using adjusted weights per cluster")
            w = (
                analysis.loc[clustering["cluster"].values, "adj_size"]
                / analysis.loc[clustering["cluster"].values, "size"]
            )
        else:
            raise ValueError("invalid value {0} for args.weights!".format(args.weights))
        w[~w.index.isin(clones) | w.index.isin(singletons)] = 0
        weights = pd.Series(w.values / w.sum(), index=np.arange(nreads))
    else:
        weights = pd.Series(np.ones(nreads) / nreads, index=np.arange(nreads))

    gap_read = gaps["read_idx"]
    gap_left = gaps["gap_left"]
    gap_right = gaps["gap_right"]
    gap_size = gaps["gap_size"]

    # breaks in consistent orientation
    Leff = Ltot // binsize
    take = (gap_size >= args.max_gap) & (gap_left // binsize < Leff) & (gap_right // binsize < Leff)
    bp_hist = scipy.sparse.csr_matrix(
        (
            weights[gap_read[take]],
            (
                np.minimum(gap_left[take], gap_right[take]) // binsize,
                np.maximum(gap_left[take], gap_right[take]) // binsize,
            ),
        ),
        shape=(Leff, Leff),
    ).todense()
    np.nan_to_num(bp_hist, copy=False)

    logger.info("creating figure and saving to {0}\n".format(args.output))
    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.markers import MarkerStyle
    matplotlib.rcParams['pdf.fonttype'] = 42

    major_ticks = []
    minor_ticks = []
    minor_labels = []
    for rec in anno_recs:
        start = shift_coord(int(rec[3][1]), cov_int)
        end = shift_coord(int(rec[3][2]), cov_int)
        major_ticks += [start, end]
        minor_ticks.append((start + end) / 2)
        minor_labels.append(rec[3][3])

    fig = plt.figure(figsize=(3, 3.75))

    #fig.text(0.535, 0.99, args.sample, size="large", ha="center", va="top")

    scale_factor = args.scale_factor
    assert Leff % scale_factor == 0, "Leff is not a multiple of scale_factor"
    bph_p = scipy.sparse.csr_matrix(
        np.asarray(bp_hist.T)
        .reshape((Leff, Leff // scale_factor, scale_factor))
        .sum(-1)
        .reshape((Leff // scale_factor, scale_factor, Leff // scale_factor))
        .sum(1)
    )
    bph_p.eliminate_zeros()

    ax = fig.add_axes([0.12, 0.03, 0.85, 0.7])
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="gray", lw=0.5, zorder=1)
    ax.hlines(np.array(major_ticks) // (binsize * scale_factor),
              #xmin=0,
              xmin=np.array(major_ticks) // (binsize * scale_factor),
              xmax=Leff // scale_factor,
              color='lightgrey',
              lw=.5,
              zorder=1)
    ax.vlines(np.array(major_ticks) // (binsize * scale_factor), 
              ymin=0,
              ymax=np.array(major_ticks) // (binsize * scale_factor),
              #ymax=Leff // scale_factor,
              color='lightgrey',
              lw=.5,
              zorder=1)
    ax.scatter(
        bph_p.nonzero()[0] + .5,
        bph_p.nonzero()[1] + .5,
        c=np.log(bph_p.data),
        cmap=plt.cm.Greys,
        marker="s",
        s=13.5,
        linewidths=0.2,
        edgecolors="k",
        zorder=2
    )
    ax.set_xlim([Leff // scale_factor, 0])
    ax.set_ylim([Leff // scale_factor, 0])
    ax.set_xticks(np.array(major_ticks)[::-1] // (binsize * scale_factor))
    ax.set_xticks(np.array(minor_ticks)[::-1] // (binsize * scale_factor), minor=True)
    ax.set_xticklabels([])
    ax.xaxis.tick_top()
    ax.set_xticklabels(minor_labels[::-1], minor=True, rotation=90)#, size='small')
    ax.set_yticks(np.array(major_ticks)[::-1] // (binsize * scale_factor))
    ax.set_yticks(np.array(minor_ticks)[::-1] // (binsize * scale_factor), minor=True)
    ax.set_yticklabels([])
    ax.set_yticklabels(minor_labels[::-1], minor=True)#, size='small')
    ax.tick_params(which="minor", length=0)
    #ax.grid(alpha=0.5, color="lightgrey", which="major")
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = fig.add_axes([0.12, 0.8, 0.85, 0.15])
    ax.plot(np.arange(Leff, (bp_hist + bp_hist.T).mean(0).A1, 'k-', lw=.5)
    #ax.plot(np.arange(Leff), bp_hist.mean(0).A1, "-", color='#900C3F',lw=0.5)
    #ax.plot(np.arange(Leff), bp_hist.T.mean(0).A1, "-", color='#E3963E',lw=0.5)
    ax.set_xlim([Leff, 0])
    ax.set_xticks(np.array(major_ticks) // (binsize))
    ax.set_xticks(np.array(minor_ticks) // (binsize), minor=True)
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.tick_params(which="minor", length=0)
    #ax.grid(alpha=0.5, color="lightgrey", which="major")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    fig.savefig(args.output, dpi=300)
