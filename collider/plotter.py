"""contains plotting functions"""
import awkward as ak
import matplotlib.pyplot as plt


def plot_histogram(variable, n_bins, range_tuple, save_name=None):
    """e.g. plot_histogram(data['lep_pt'], 100, (0, 100))"""
    all_lep_pt = ak.to_numpy(ak.flatten(variable, axis=None))
    plt.hist(all_lep_pt, n_bins, range_tuple)
    if save_name is None:
        plt.show()
    else:
        plt.savefig(save_name)
        plt.close()
