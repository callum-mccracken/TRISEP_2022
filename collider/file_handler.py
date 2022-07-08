"""file_handler.py handles reading data from root files"""

import os
from typing import List
import uproot

OPEN_DATA_URL = 'http://opendata.atlas.cern/release/samples/Data/'  # +filename


def open_file(filename: str,
              branches:List[str]=None,
              sample_size:int=None):
    """
    Open a root file, for the given particle type, branches, and sample size

    filename: something like 'DataMuons.root'
    branches: branches you need: ['lep_pt', 'lep_eta', ...], less = quicker
    sample_size: how many events to return, use smaller samples for testing

    """
    # get filename, either as the file in this directory or the OpenData link
    if not os.path.exists(filename):
        filename = OPEN_DATA_URL + filename

    # open the file
    root_file = uproot.open(filename)

    # Load the 'mini' tree within the file
    tree = root_file["mini"]

    # get branches we care about, as AwkwardArrays
    # specifying branches makes this faster, not technically needed
    if branches is None:
        arrays = tree.arrays()
    else:
        arrays = tree.arrays(branches)

    # Use a smaller dataset for testing purposes!
    if sample_size is None:
        return arrays
    return arrays[:sample_size]


if __name__ == '__main__':
    open_file('DataMuons.root')  # takes like 10 seconds to run
