# Classes for handling the data 
import numpy as np
import pandas as pd
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(os.path.dirname(FILE_PATH), "data")

class TurbulentFlatData:
    """Base class for handling turbulent flat channel data."""
    def __init__(self, re_number):
        self.re_number = re_number
        self.mean_velocity = None
        self.rms_velocities = {}
    
    def load_mean_velocity(self):
        """Load mean velocity data."""
        filepath = os.path.join(DATA_PATH, f"Re{self.re_number}_MeanVelocity.csv")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        self.mean_velocity = pd.read_csv(filepath)
        return self.mean_velocity
    
    def load_rms_velocities(self):
        """Load RMS velocity data for all components."""
        for component in ['U', 'V', 'W']:
            filepath = os.path.join(DATA_PATH, f"Re{self.re_number}_RMS_{component}.csv")
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"File not found: {filepath}")
            self.rms_velocities[component] = pd.read_csv(filepath)
        return self.rms_velocities
    
    def get_mean_velocity(self):
        """Get mean velocity data, loading if necessary."""
        if self.mean_velocity is None:
            self.load_mean_velocity()
        return self.mean_velocity
    
    def get_rms_velocity(self, component):
        """Get RMS velocity data for specified component."""
        if component not in self.rms_velocities:
            if not self.rms_velocities:  # If empty, load all
                self.load_rms_velocities()
        return self.rms_velocities.get(component)

class re_180(TurbulentFlatData):
    def __init__(self):
        super().__init__(180)