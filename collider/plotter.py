"""contains plotting functions"""
from typing import Tuple
import awkward as ak
import matplotlib.pyplot as plt


def plot_histogram(variable: ak.Array,
                   n_bins: int,
                   range_tuple: Tuple[float],
                   save_name: str=None,
                   xlabel: str=None,
                   ylabel: str=None,
                   title: str=None):
    """e.g. plot_histogram(data['lep_pt'], 100, (0, 100))"""
    all_lep_pt = ak.to_numpy(ak.flatten(variable, axis=None))
    plt.hist(all_lep_pt, n_bins, range_tuple)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)    

    if save_name is None:
        plt.show()
    else:
        plt.savefig(save_name)
        plt.close()
