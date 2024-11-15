import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser=ArgumentParser()

parser.add_argument("--figure", dest="figure", help="""output figure""")
parser.add_argument("--samples", dest="samples", help="""list of samples""")
parser.add_argument("--color", dest="color", default='cluster', help="""what to color circles by [cluster]""")
parser.add_argument(
    "--fig_width",
    dest="fig_width",
    default=3,
    type=float,
    help="""fig width""",
)
parser.add_argument(
    "--fig_height",
    dest="fig_height",
    type=float,
    default=.75,
    help="""figure height in inches""",
)
parser.add_argument(
    "--dpi",
    dest="dpi",
    default=300,
    type=int,
    help="""dpi of output figure [300]""",
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

args = parser.parse_args()

samples = args.samples.split(',')

figsize = (args.fig_width, args.fig_height)
fig, axs = plt.subplots(1, len(samples), figsize=figsize, squeeze=True, sharey=False)
fig.subplots_adjust(bottom=.5)#, left=.2)

isotype_colors = {"SM":"#000000", "SG3":"#98E4FD", "SG1":"#25C5C5","SA1":"#D92916", 
                  "SG2":"#147F7F", "SG4":"#084C4C", "SE":"#C9C9C9", "SA2":"#8A170C"}

for ns, sample in enumerate(samples):
    logger.info("loading results for {0}".format(sample))
    clustering = pd.read_csv(os.path.join('pipeline',sample,sample + '_clustering.csv'), 
                             index_col=0, header=0)
    cluster_analysis = pd.read_csv(os.path.join('pipeline',sample,sample+'_cluster_analysis.csv'), 
                                   index_col=0, header=0)

    clusters = clustering["filtered_cluster"].dropna()
    clones = clusters[clusters >= 0].astype(int).unique()
    cluster_isotype_count = cluster_analysis.loc[clones, "isotype"].dropna().value_counts()
    nclusters = cluster_isotype_count.sum()

    isotype_fracs = pd.DataFrame({'1': 100 * cluster_isotype_count / nclusters})

    isotype_fracs.T.plot(kind="barh", stacked=True, ax=axs[ns],
                         color=[isotype_colors[it] for it in isotype_fracs.index], legend=False)
    
    axs[ns].set_xlabel("% clusters", size='small')
    axs[ns].set_xticks([0,50,100])
    axs[ns].set_xticklabels(['0','50','100'], size='small')
    axs[ns].set_yticks([])
    axs[ns].set_yticklabels([])
    axs[ns].spines['right'].set_visible(False)
    axs[ns].spines['top'].set_visible(False)
    axs[ns].spines['left'].set_visible(False)
    axs[ns].set_ylim([-.3,.3])
    #if ns==len(samples) -1:
    #    axs[ns].legend(loc=2,bbox_to_anchor=(1,1), frameon=False)

plt.tight_layout()
fig.savefig(args.figure)
