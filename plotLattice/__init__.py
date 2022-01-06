import numpy as np
#  import matplotlib.pyplot as plt
#  import csv
#  from scipy.signal import find_peaks

# reduced Planck unit
m_pl = 2.44e18

#TODO: resolve some errors (mostly to add argument in the function)

class GlobVar:
    '''
    Class for 'global' variables, which are lattice parameters and etc
    '''
    def __init__(self, fStar, omegaStar, L, N, phi0, m_eff_min):
        self.fStar = fStar
        self.omegaStar = omegaStar
        self.L = L
        self.N = N
        self.phi0 = phi0
        self.m_eff_min = m_eff_min

    def print(self):
        print(f"Parameters are (all in m_pl):\n\tf* = {self.fStar}\n\tomega* = {self.omegaStar}\n\tL = {self.L}\n\tN = {self.N}\n\tphi0 = {self.phi0}\n\tm_eff_min = {self.m_eff_min}")


def initialize(fn):
    """
    Read parameter files, src/model/parameter-files/'fn'.
    Only need these parameters now.
    """
    with open(fn) as f:
        for line in f:
            line = line.partition('#')[0]
            line = line.rstrip()
            if line.startswith('initial_amplitudes = '):
                init_amp = float(line.replace("initial_amplitudes = ", ''))
            elif line.startswith('d = '):
                d = float(line.replace("d = ", ''))
            elif line.startswith('kIR = '):
                kIR = float(line.replace('kIR = ', ''))
                L = 2*np.pi/kIR
            elif line.startswith('lSide = '):
                L = float(line.replace('lSide = ', ''))
            elif line.startswith('N = '):
                N = int(line.replace('N = ', ''))
            elif line.startswith('kCutOff = '):
                kCutOff = float(line.replace('kCutOff = ', ''))
            elif line.startswith('phi0 = '):
                phi0 = float(line.replace('phi0 = ', ''))
    # to convert all in reduced Planck
    fStar = init_amp/m_pl
    omegaStar = np.sqrt(d)/3*phi0/m_pl
    L /= omegaStar  # originally in prog. unit
    kCutOff *= omegaStar  # originally in prog. unit
    phi0 /= m_pl  # originally in GeV
    m_eff_min = 2.97e-8 * phi0**2

    delta_x = L/N  # lattice spacing
    k_UV = np.sqrt(3)*np.pi/delta_x
    k_IR = 2*np.pi/L
    print("Momentum cutoffs are: k_IR = %e m_pl, k_UV = %e m_pl"
          % (k_IR, k_UV))
    print("Initializing momenta (UV) cutoff is: kCutoff = %e m_pl"
          % (kCutOff))
    print("")

    return GlobVar(fStar, omegaStar, L, N, phi0, m_eff_min)
