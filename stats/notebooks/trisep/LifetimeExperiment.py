# A LifetimeExperiment creates a single unstable isotope, which remains in the detector for
# at most 0.2 seconds (le_max_time). Provided it decays, the time it takes to decay is reported.
# The timing measurement is characterized by a resolution and bias.
# There is also a background process that produces uniform times (-0.2,0.2)
# Reported times are between -0.2 and 0.2 seconds (collecting most of the decays). A fixed number of
# events are returned: counting_time/le_max_time

import pickle
import time

import numpy as np
from scipy import stats


class LifetimeExperiment:
    def __init__(self, bypass_wait=False):
        self.counting_time = 10.
        self.bypass_wait = bypass_wait
        with open('trisep.p', 'rb') as f:
            TC = pickle.load(f)
        self.le_max_time = 0.2
        self.le_isotope_lifetime = TC.le_isotope_lifetime  # should be much less than le_max_time
        self.le_time_resolution = TC.le_time_resolution
        self.le_time_offset = TC.le_time_offset
        self.le_background_fraction = TC.le_background_fraction  # fraction of data coming from background
        self.times = None

        print("Lab lifetime experiment built. Default counting time is", self.counting_time, "seconds.")

    def set_counting_time(self, counting_time):
        if counting_time > 30.:
            print('Counting time not changed:', counting_time, 'seconds is too long to wait! (30 seconds max)')
        elif counting_time < 0.:
            print('Error: Counting time must not be negative')
        else:
            self.counting_time = counting_time

    def get_counting_time(self):
        return self.counting_time

    def get_times(self):
        return self.times

    def start(self):
        if not self.bypass_wait:
            print('Please wait', self.counting_time, 'seconds...')
            time.sleep(self.counting_time)
        self.produce_times()

    def produce_times(self):
        n_events = int(self.counting_time / self.le_max_time)
        times = []
        for i in range(n_events):
            if stats.uniform.rvs() < self.le_background_fraction:
                # a background event
                uniform_time = self.le_max_time * 2. * (stats.uniform.rvs() - 0.5)
                times.append(uniform_time)
            else:
                # produce a real event within the acceptable range
                in_range = False
                observed_time = None
                while not in_range:
                    decay_time = stats.expon.rvs(scale=self.le_isotope_lifetime)
                    # add offset and resolution
                    observed_time = decay_time + stats.norm.rvs(self.le_time_offset, self.le_time_resolution)
                    in_range = -self.le_max_time < observed_time < self.le_max_time
                times.append(observed_time)
        self.times = times
        return


class SimulatedLifetimeExperiment(LifetimeExperiment):
    def __init__(self, isotope_lifetime=0.05, time_resolution=0., time_offset=0., background_fraction=0.):
        self.counting_time = 10.
        self.le_max_time = 0.2
        self.le_isotope_lifetime = isotope_lifetime
        self.le_time_resolution = time_resolution
        self.le_time_offset = time_offset
        self.le_background_fraction = background_fraction
        print("Simulated lifetime experiment built. Counting time =", self.counting_time,
              "Time resolution =", time_resolution, "Time offset=", time_offset,
              "Background fraction =", background_fraction)
        self.times = None

    def set_counting_time(self, counting_time):
        if counting_time < 0.:
            print('Error: Counting time must not be negative')
        else:
            self.counting_time = counting_time

    def set_isotope_lifetime(self, isotope_lifetime):
        if 0. <= isotope_lifetime <= self.le_max_time:
            self.le_isotope_lifetime = isotope_lifetime
        else:
            print('Error: Isotope lifetime must be between 0. and', self.le_max_time)

    def set_time_resolution(self, time_resolution):
        if 0. <= time_resolution <= self.le_max_time:
            self.le_time_resolution = time_resolution
        else:
            print('Error: Time resolution must be between 0. and', self.le_max_time)

    def set_time_offset(self, time_offset):
        if -self.le_max_time <= time_offset <= self.le_max_time:
            self.le_time_offset = time_offset
        else:
            print('Error: Time offset must be between', -self.le_max_time, 'and', self.le_max_time)

    def set_background_fraction(self, background_fraction):
        if 0. <= background_fraction <= 1.:
            self.le_background_fraction = background_fraction
        else:
            print('Error: Background fraction must be between 0. and 1.')

    def start(self):
        self.produce_times()

    def get_pdf(self, observed_time):
        # Calculate the pdf
        bf = self.le_background_fraction
        uniform_pdf = 0.
        special_pdf = 0.
        if -self.le_max_time <= observed_time < self.le_max_time:
            # Uniform
            uniform_pdf = 1. / (2. * self.le_max_time)

            # Convolution of exponential and normal (CUTOFF AT +/- le_max_time, so MUST renormalize)
            # catch issues where pdf might be very small
            if observed_time - self.le_time_offset > -5. * self.le_time_resolution:
                if self.le_time_resolution == 0.:  # pure exponential
                    # integral over all possible observations:
                    cdfs = stats.expon.cdf([-self.le_max_time, self.le_max_time], self.le_time_offset,
                                           self.le_isotope_lifetime)
                    integral = cdfs[1] - cdfs[0]

                    special_pdf = stats.expon.pdf(observed_time, self.le_time_offset,
                                                  self.le_isotope_lifetime) / integral
                else:  # exponential convoluted with normal
                    # integral over possible range
                    cdfs = stats.exponnorm.cdf([-self.le_max_time, self.le_max_time],
                                               self.le_isotope_lifetime / self.le_time_resolution,
                                               self.le_time_offset, self.le_time_resolution)
                    integral = cdfs[1] - cdfs[0]

                    special_pdf = stats.exponnorm.pdf(observed_time, self.le_isotope_lifetime / self.le_time_resolution,
                                                      self.le_time_offset, self.le_time_resolution) / integral

        return bf * uniform_pdf + (1. - bf) * special_pdf

    def get_log_likelihood(self, times):
        log_lik = -np.Inf
        zero_pdf = False

        sum_log_pdf = 0.
        for observed_time in times:
            pdf = self.get_pdf(observed_time)
            if pdf <= 0.:
                zero_pdf = True
                break
            sum_log_pdf += np.log(pdf)

        if not zero_pdf:
            log_lik = sum_log_pdf

        return log_lik
