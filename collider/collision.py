"""for doing the main thing -- TODO fix this docstring"""

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
plot_histogram(data['lep_pt'] / 1000, 25, (0, 100), "lep_pt_without_cuts.png")
plot_histogram(data['lep_n'], 3, (1, 4), "lep_n_before_cuts.png")
# 2nd one justifies cutting on exactly two -- only lose 85 events

# make cuts
data = cuts.make_cuts(data, particle_id=PARTICLE_ID)
print(len(data), 'events left after cuts')

# get data in pretty variables
lep_pt = data['lep_pt'] / 1000  # to convert MeV to GeV
lep_eta = data['lep_eta']
lep_phi = data['lep_phi']
lep_E = data['lep_E'] / 1000  # MeV --> GeV
lep_n = data['lep_n']

# plot combined lepton properties after cuts
plot_histogram(lep_pt, 25, (0, 100), "lep_pt.png")
plot_histogram(lep_eta, 25, (0, 4), "lep_eta.png")
plot_histogram(lep_phi, 25, (0, np.pi), "lep_phi.png")
plot_histogram(lep_E, 25, (0, 100), "lep_E.png")

# and individual properties
plot_histogram(lep_pt[:, 0], 25, (0, 100), "lep_0_pt.png")
plot_histogram(lep_eta[:, 0], 25, (0, 4), "lep_0_eta.png")
plot_histogram(lep_phi[:, 0], 25, (0, np.pi), "lep_0_phi.png")
plot_histogram(lep_E[:, 0], 25, (0, 100), "lep_0_E.png")

plot_histogram(lep_pt[:, 1], 25, (0, 100), "lep_1_pt.png")
plot_histogram(lep_eta[:, 1], 25, (0, 4), "lep_1_eta.png")
plot_histogram(lep_phi[:, 1], 25, (0, np.pi), "lep_1_phi.png")
plot_histogram(lep_E[:, 1], 25, (0, 100), "lep_1_E.png")



# make vectors, construct Z candidates
lepton_vectors = vector.zip({
    'pt': lep_pt, 'eta': lep_eta, 'phi': lep_phi, 'E': lep_E})
lead_lepton = lepton_vectors[:, 0]
next_lepton = lepton_vectors[:, 1]
z_candidate = lead_lepton + next_lepton

# reconstructed mass plot
plot_histogram(z_candidate.mass, 100, (0, 100), "reco_z_mass.png")


plot_histogram(z_candidate.pt, 100, (0, 100), "reco_z_pt.png")
plot_histogram(z_candidate.eta, 100, (0, 100), "reco_z_eta.png")
plot_histogram(z_candidate.phi, 100, (0, 100), "reco_z_phi.png")
plot_histogram(z_candidate.E, 100, (0, 100), "reco_z_E.png")
