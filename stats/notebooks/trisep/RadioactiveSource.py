# A RadioactiveSource contains unstable isotopes
import pickle

from scipy import stats


class RadioactiveSource:
    def __init__(self, activity):
        self.activity = None
        with open('trisep.p', 'rb') as f:
            TC = pickle.load(f)
        self.max_calibration_source_activity = TC.max_calibration_source_activity
        self.set_activity(activity)

    def set_activity(self, activity):
        if activity < 0.:
            print('Error: source activity must not be negative!')
        else:
            self.activity = activity
            if activity > self.max_calibration_source_activity:
                print('Warning: the activity specified is higher than the lab rules allow (',
                      self.max_calibration_source_activity, 'Bq)')

    def get_decays(self, recording_time):
        number = 0
        if recording_time >= 0.:
            expected_value = self.activity * recording_time
            number = stats.poisson.rvs(expected_value)
        else:
            print('Error: recording time must not be negative!')
        return number
