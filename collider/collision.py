"""Main module for the collider physics tutorial at TRISEP 2022"""

# physics
import numpy as np
import vector
import file_handler
import cuts
from plotter import plot_histogram

# constants and such
BRANCH_NAMES = [  # things we'll need from file (use None to get everything)
    'lep_pt',
    'lep_eta',
    'lep_phi',
    'lep_E',
    'lep_n',
    'lep_charge',
    'lep_type']
SAMPLE_SIZE = 100000  # use None to get the whole dataset

FILE_NAME = 'DataMuons.root'
PARTICLE_ID = {
    # "file_name": PDG ID,
    "DataMuons.root": 13,
    "DataEgamma.root": 11,
}[FILE_NAME]

# open file, get data
data = file_handler.open_file(
    FILE_NAME, branches=BRANCH_NAMES, sample_size=SAMPLE_SIZE)
print(len(data), 'events to start')

# make some pre-cut histograms
plot_histogram(
    data['lep_pt'] / 1000, n_bins=25, range_tuple=(0, 100),
    xlabel="pT [GeV]", ylabel="counts", title="Lepton pT before cuts",
    save_name="lep_pt_without_cuts.png")
plot_histogram(
    data['lep_n'], n_bins=3, range_tuple=(1, 4),
    xlabel="lepton number", ylabel="counts", title="Lepton number before cuts",
    save_name="lep_n_before_cuts.png")
# 2nd one justifies cutting on exactly two -- only lose 85 events

# make cuts (2 leptons, same type opposite sign, pt>25GeV)
data = cuts.make_cuts(data, particle_id=PARTICLE_ID)
print(len(data), 'events left after cuts')

# get data in pretty variables
lep_pt = data['lep_pt'] / 1000  # to convert MeV to GeV
lep_eta = data['lep_eta']
lep_phi = data['lep_phi']
lep_E = data['lep_E'] / 1000  # MeV --> GeV
lep_n = data['lep_n']

# plot combined lepton properties after cuts
plot_histogram(
    lep_pt, n_bins=25, range_tuple=(0, 100),
    xlabel="lepton pT", ylabel="counts", title="Lepton pT after cuts",
    save_name="lep_pt.png")
plot_histogram(
    lep_eta, n_bins=25, range_tuple=(0, 4),
    xlabel="lepton eta", ylabel="counts", title="Lepton eta after cuts",
    save_name="lep_eta.png")
plot_histogram(
    lep_phi, n_bins=25, range_tuple=(0, np.pi),
    xlabel="lepton phi", ylabel="counts", title="Lepton phi after cuts",
    save_name="lep_phi.png")
plot_histogram(
    lep_E, n_bins=25, range_tuple=(0, 100),
    xlabel="lepton E", ylabel="counts", title="Lepton E after cuts",
    save_name="lep_E.png")

# and individual properties
plot_histogram(
    lep_pt[:, 0], n_bins=25, range_tuple=(0, 100),
    xlabel="leading lepton pT", ylabel="counts",
    title="Leading lepton pT after cuts",
    save_name="lep_0_pt.png")
plot_histogram(
    lep_eta[:, 0], n_bins=25, range_tuple=(0, 4),
    xlabel="leading lepton eta", ylabel="counts",
    title="Leading lepton eta after cuts",
    save_name="lep_0_eta.png")
plot_histogram(
    lep_phi[:, 0], n_bins=25, range_tuple=(0, np.pi),
    xlabel="leading lepton phi", ylabel="counts",
    title="Leading lepton phi after cuts",
    save_name="lep_0_phi.png")
plot_histogram(
    lep_E[:, 0], n_bins=25, range_tuple=(0, 100),
    xlabel="leading lepton E", ylabel="counts",
    title="Leading lepton E after cuts",
    save_name="lep_0_E.png")

plot_histogram(
    lep_pt[:, 1], n_bins=25, range_tuple=(0, 100),
    xlabel="sub-leading lepton pT", ylabel="counts",
    title="sub-leading lepton pT after cuts",
    save_name="lep_1_pt.png")
plot_histogram(
    lep_eta[:, 1], n_bins=25, range_tuple=(0, 4),
    xlabel="sub-leading lepton eta", ylabel="counts",
    title="sub-leading lepton eta after cuts",
    save_name="lep_1_eta.png")
plot_histogram(
    lep_phi[:, 1], n_bins=25, range_tuple=(0, np.pi),
    xlabel="sub-leading lepton phi", ylabel="counts",
    title="sub-leading lepton phi after cuts",
    save_name="lep_1_phi.png")
plot_histogram(
    lep_E[:, 1], n_bins=25, range_tuple=(0, 100),
    xlabel="sub-leading lepton E", ylabel="counts",
    title="sub-leading lepton E after cuts",
    save_name="lep_1_E.png")

# make vectors, construct Z candidates
lepton_vectors = vector.zip({
    'pt': lep_pt, 'eta': lep_eta, 'phi': lep_phi, 'E': lep_E})
lead_lepton = lepton_vectors[:, 0]
next_lepton = lepton_vectors[:, 1]
z_candidate = lead_lepton + next_lepton

# plot reconstructed Z properties
plot_histogram(
    z_candidate.mass, n_bins=100, range_tuple=(50, 150),
    xlabel="reconstructed Z mass", ylabel="counts",
    title="reconstructed Z mass",
    save_name="reco_z_mass.png")
plot_histogram(
    z_candidate.pt, n_bins=25, range_tuple=(0, 100),
    xlabel="reconstructed Z pt", ylabel="counts",
    title="reconstructed Z pt",
    save_name="reco_z_pt.png")
plot_histogram(
    z_candidate.eta, n_bins=25, range_tuple=(0, 4),
    xlabel="reconstructed Z eta", ylabel="counts",
    title="reconstructed Z eta",
    save_name="reco_z_eta.png")
plot_histogram(
    z_candidate.phi, n_bins=25, range_tuple=(0, np.pi),
    xlabel="reconstructed Z phi", ylabel="counts",
    title="reconstructed Z phi",
    save_name="reco_z_phi.png")
plot_histogram(
    z_candidate.E, n_bins=25, range_tuple=(0, 100),
    xlabel="reconstructed Z E", ylabel="counts",
    title="reconstructed Z E",
    save_name="reco_z_E.png")
