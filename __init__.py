import numpy as np
import matplotlib.pyplot as plt
import csv
#  from scipy.signal import find_peaks

plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# reduced Planck unit
m_pl = 2.44e18


def read_para(fn):
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
    # all in reduced Planck
    fStar = init_amp/m_pl
    omegaStar = np.sqrt(d)/3*phi0/m_pl
    L /= omegaStar
    kCutOff *= omegaStar
    phi0 /= m_pl
    parameters = {
        "f*": fStar,
        "omega*": omegaStar,
        "L": L,
        "N": N,
        "kCutOff": kCutOff,
        "phi0": phi0
    }
    print("Lattice parameters (all in m_pl):")
    print(parameters)
    return parameters


#  def get_hubble_inf(phi0):
    #  '''
    #  return the scale of inflation in SFPI model fitted with Planck
    #  '''
    #  H = 8.6e-9 * phi0**3  # m_pl
    #  return H


def initialize(parameter_fn):
    """
    Read lattice parameters from parameter_fn and save global variables
    """
    para_dict = read_para("./sfpi.in")
    global fStar, omegaStar, L, N, phi0, H0, m_eff_min
    fStar = para_dict["f*"]
    omegaStar = para_dict["omega*"]
    L = para_dict["L"]
    N = para_dict['N']
    phi0 = para_dict['phi0']
    #  H0 = get_hubble_inf(phi0)
    m_eff_min = 2.97e-8 * phi0**2  # in m_pl

    delta_x = L/N
    k_UV = np.sqrt(3)*np.pi/delta_x
    k_IR = 2*np.pi/L
    print("Momentum cutoffs are: k_IR = %e m_pl, k_UV = %e m_pl"
          % (k_IR, k_UV))
    kCutOff = para_dict["kCutOff"]
    print("Initializing momenta (UV) cutoff is: kCutoff = %e m_pl"
          % (kCutOff))
    print("")
