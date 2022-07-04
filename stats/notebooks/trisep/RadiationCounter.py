# A RadiationCounter accepts a RadiationSource and counts the number of decays over a fixed period of time
import pickle
import time

from scipy import stats

from trisep.RadioactiveSource import RadioactiveSource


class RadiationCounter:
    def __init__(self, bypass_wait=False):
        self.counting_time = 10.
        self.bypass_wait = bypass_wait
        with open('trisep.p', 'rb') as f:
            TC = pickle.load(f)
        self.max_calibration_source_activity = TC.max_calibration_source_activity
        self.lab_source = RadioactiveSource(TC.source_activity)
        self.source = None
        self.efficiency = TC.counter_efficiency
        self.background = TC.counter_background
        self.count = 0
        print("Lab radiation counter built. Default counting time is", self.counting_time, "seconds.")

    def set_counting_time(self, counting_time):
        if counting_time > 30.:
            print('Counting time not changed:', counting_time, 'seconds is too long to wait! (30 seconds max)')
        elif counting_time < 0.:
            print('Error: Counting time must not be negative')
        else:
            self.counting_time = counting_time

    def get_counting_time(self):
        return self.counting_time

    def insert_lab_source(self):
        self.source = self.lab_source

    def insert_known_source(self, source):
        if source.activity > self.max_calibration_source_activity:
            print('Error: The known source activity is higher than lab regulations')
        else:
            self.source = source

    def remove_source(self):
        self.source = None

    def get_count(self):
        return self.count

    def start(self):
        if not self.bypass_wait:
            print('Please wait', self.counting_time, 'seconds...')
            time.sleep(self.counting_time)
        self.produce_count(verbose=True)

    def produce_count(self, verbose=False):
        # A single Poisson random variable would be faster but
        # the following is written to be more transparent

        # Number of decays during the recording time:
        n_decays = 0
        if self.source is not None:
            n_decays = self.source.get_decays(self.counting_time)
        elif verbose:
            print('Warning: There is no source in the detector')

        # Number of signal events: accounting for efficiency
        n_signal = n_decays
        if self.efficiency != 1.:
            n_signal = stats.binom.rvs(n_decays, self.efficiency)

        # Number of background events during the recording time:
        expected_background = self.background * self.counting_time
        n_background = stats.poisson.rvs(expected_background)

        self.count = n_signal + n_background
        return


class SimulatedRadiationCounter(RadiationCounter):
    def __init__(self, efficiency=1., background=0.):
        self.counting_time = 10.
        self.source = None
        self.efficiency = efficiency
        self.background = background
        print("Simulated detector built. Counting time =", self.counting_time,
              "Efficiency =", efficiency, "Background rate=", background, "(Hz)")

    def set_counting_time(self, counting_time):
        if counting_time < 0.:
            print('Error: Counting time must not be negative')
        else:
            self.counting_time = counting_time

    def insert_lab_source(self):
        print('Error: You cannot put the lab source into a simulated detector!')

    def insert_known_source(self, source):
        self.source = source

    def set_efficiency(self, efficiency):
        if 0. <= efficiency <= 1.:
            self.efficiency = efficiency
        else:
            print('Error: Efficiency must be between 0. and 1.')

    def set_background(self, background):
        if background >= 0.:
            self.background = background
        else:
            print('Error: Background rate cannot be negative')

    def start(self):
        self.produce_count(verbose=False)

    def get_likelihood(self, count):
        activity = 0
        if self.source is not None:
            activity = self.source.activity
        expected_value = activity * self.counting_time * self.efficiency + self.background * self.counting_time
        likelihood = stats.poisson.pmf(count, expected_value)

        return likelihood
