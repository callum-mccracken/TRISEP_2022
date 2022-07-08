"""a module to store the cuts function"""
import awkward

def make_cuts(data: awkward.Array, particle_id: int):
    """
    A function to make cuts based on what particles we have.

    data: unfiltered data from one of the files,
    particle_id: PDG ID of the lepton used in the file, 11/13 for e/mu"""
    # Ensure we have at least two leptons (exactly 2)
    data = data[data['lep_n'] == 2]

    # Ensure the leptons are of the right type
    data = data[(data['lep_type'][:, 0] == particle_id) *
                (data['lep_type'][:, 0] == particle_id)]

    # Ensure the leptons have the opposite charge (since they came from a Z)
    data = data[data['lep_charge'][:, 0] == -data['lep_charge'][:, 1]]

    # Ensure both leptons have pT >= 25GeV
    data = data[(data['lep_pt'][:, 0] >= 25000) *
                (data['lep_pt'][:, 1] >= 25000)]

    # more cuts here if needed ...
    return data
