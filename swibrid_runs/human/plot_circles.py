import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser=ArgumentParser()

parser.add_argument("--figure", dest="figure", help="""output figure""")
parser.add_argument("--samples", dest="samples", help="""list of samples""")
parser.add_argument("--color", dest="color", default='cluster', help="""what to color circles by [cluster]""")
parser.add_argument("--ncols", dest="ncols", default=3, type=int, help="""number of columns [3]""")
parser.add_argument(
    "--fig_width",
    dest="fig_width",
    default=1,
    type=float,
    help="""width of panels in inches [5]""",
)
parser.add_argument(
    "--fig_height",
    dest="fig_height",
    type=float,
    default=1,
    help="""width of panels in inches [4.5]""",
)
parser.add_argument(
    "--dpi",
    dest="dpi",
    default=300,
    type=int,
    help="""dpi of output figure [300]""",
)
parser.add_argument(
    "--use_circlify",
    dest="use_circlify",
    default=False,
    action="store_true",
    help="""use circify package for circular packing""",
)

def rand_cmap(
    nlabels,
    type="bright",
    first_color_black=False,
    last_color_black=False,
    verbose=False,
    bad="gray",
    under="gray",
    over="gray",
    seed=0,
):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    (from https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib/14720445)
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    np.random.seed(seed)

    # Generate color map for bright colors, based on hsv
    if type == "bright":
        randHSVcolors = [
            (
                np.random.uniform(low=0.0, high=1),
                np.random.uniform(low=0.2, high=1),
                np.random.uniform(low=0.9, high=1),
            )
            for i in range(nlabels)
        ]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

    # Generate soft pastel colors, by limiting the RGB spectrum
    elif type == "soft":
        low = 0.6
        high = 0.95
        randRGBcolors = [
            (
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
                np.random.uniform(low=low, high=high),
            )
            for i in range(nlabels)
        ]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

    else:
        print('Please choose "bright" or "soft" for type')
        return

    random_colormap = LinearSegmentedColormap.from_list("new_map", randRGBcolors, N=nlabels)

    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        colorbar.ColorbarBase(
            ax,
            cmap=random_colormap,
            norm=norm,
            spacing="proportional",
            ticks=None,
            boundaries=bounds,
            format="%1i",
            orientation="horizontal",
        )

    random_colormap.set_bad(bad)
    random_colormap.set_under(under)
    random_colormap.set_over(over)

    return random_colormap


import os
import sys
import numpy as np
import re
import pandas as pd
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from logzero import logger

matplotlib.rcParams.update({"font.size": 8})
matplotlib.rcParams.update({"axes.linewidth": 0.5})
matplotlib.rcParams.update({"xtick.major.width": 0.5})
matplotlib.rcParams.update({"ytick.major.width": 0.5})
matplotlib.rcParams['pdf.fonttype'] = 42

isotype_colors = {"SM":"#000000", "SG3":"#98E4FD", "SG1":"#25C5C5","SA1":"#D92916",
                  "SG2":"#147F7F", "SG4":"#084C4C", "SE":"#C9C9C9", "SA2":"#8A170C"}

args = parser.parse_args()

samples = args.samples.split(',')
figsize = (args.fig_width * args.ncols, args.fig_height * np.ceil(len(samples) / args.ncols))
fig, axs = plt.subplots(int(np.ceil(len(samples)/args.ncols)), args.ncols, 
                        figsize=figsize, sharex=True, sharey=True, squeeze=False)

fig.subplots_adjust(top=.99,bottom=.01, left=.01, right=.99, hspace=0, wspace=0)

lims=[]
for ns, sample in enumerate(samples):
    clustering_results = os.path.join('pipeline',sample,sample + '_clustering.csv')
    logger.info("loading clustering from {0}".format(clustering_results))
    clustering = pd.read_csv(clustering_results, index_col=0, header=0)
    summary = pd.read_csv(os.path.join('pipeline',sample,sample+'_summary.csv'), index_col=0).squeeze()
    reads = clustering["cluster"].dropna().index
    nreads = len(reads)
    clustering = clustering.loc[reads]
    csize = clustering.groupby(['cluster', 'filtered_cluster']).size().reset_index(level=1).sort_values([0, 'filtered_cluster'], ascending=False).drop('filtered_cluster', axis=1).squeeze()
    nclust = len(csize)
    if args.color == 'cluster':
        cmap = rand_cmap(
            max(2, nclust),
            type="bright",
            first_color_black=False,
            last_color_black=False,
            verbose=False,
            seed=10,
        )
        cluster_color = dict((k,cmap(n)) for n,k in enumerate(csize.index))
    else:
        cluster_color = clustering[['cluster',args.color]].drop_duplicates().set_index('cluster').squeeze().astype('category')
        if args.color == 'isotype':
            cluster_color = dict((n,isotype_colors[k]) for n,k in cluster_color.items())
        else:
            cluster_color = dict((n,plt.cm.Set1(k)) for n,k in cluster_color.cat.codes.items())

    cluster_alpha = clustering[['cluster','filtered_cluster']].drop_duplicates().set_index('cluster').squeeze()
    cluster_alpha = dict((n,1 if k >= 0 else 0.25) for n,k in cluster_alpha.items())

    stat_string = "{0}: {1} reads\nnclusters: {2:.0f}; nclusters_final: {3:.0f}; nclusters_eff: {4:.2f}\nentropy: {5:.2f}, inverse simpson: {6:.2f}"
    stats = stat_string.format(sample, nreads, summary['nclusters_initial'],
                               summary['nclusters_final'], summary['nclusters_eff'],
                               summary['cluster_entropy'], summary['cluster_inverse_simpson'])

    if nclust < 500 and args.use_circlify:
        import circlify

        circles = dict(zip(csize.index,
                           circlify.circlify(
                    csize.tolist(),
                    show_enclosure=False,
                    target_enclosure=circlify.Circle(x=0, y=0, r=1),
                    )))
    else:
        import packcircles

        circles = dict(zip(csize.index[::-1],list(packcircles.pack(np.sqrt(csize)[::-1].tolist()))))

    ax=axs[ns//args.ncols, ns%args.ncols]
    ax.axis("off")

    lim = 0
    for n in csize.index:
        x, y, r = circles[n]
        ax.add_patch(
            plt.Circle(
                (x, y),
                r,
                facecolor=cluster_color[n],
                edgecolor="k",
                alpha=cluster_alpha[n],
                linewidth=0.1,
                clip_on=False,
            )
        )
        m = max(abs(x) + r, abs(y) + r)
        if m > lim:
            lim = m

    lims.append(lim)
    #ax.set_title(stats, size='small', y=.95)

lim = np.max(lims)
for ns, sample in enumerate(samples):
    ax=axs[ns//args.ncols, ns%args.ncols]
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_aspect('equal')

for ns in range(len(samples), np.prod(axs.shape)):
    axs[ns//args.ncols, ns%args.ncols].axis("off")

fig.savefig(args.figure)
