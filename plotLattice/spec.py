import numpy as np
import matplotlib.pyplot as plt
import csv

plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# import global variables
#  from . import fStar, omegaStar, L, N, phi0, m_eff_min


def convert_line_2_floats(s):
    """
    Interpret line of floats separated by space as floats
    """
    floats = [float(value) for value in s.split(' ')]
    return floats


def parse_spec(file_name):
    """
    Read spectra*.dat file from CosmoLattice
    Format is for each line:
        k, P_phi, P_phi', n, Delta n_bin
    Data at different times are separated by a blanck line.
    Which time is in spectratimes*.dat
    Return:
        Big array, 0-axis is time, 1- and 2-axis are different k and spectra
    """
    spectra = None
    temp_array = None
    with open(file_name, mode='r') as file:
        csv_file = csv.reader(file)
        for line in csv_file:
            if line != []:  # if data
                if temp_array is None:  # create temp_array if non-existent
                    temp_array = np.array([convert_line_2_floats(line[0])])
                else:
                    temp_array = np.append(temp_array, np.array([convert_line_2_floats(line[0])]), axis=0)
                # keep program units
                """
                temp_array[-1, 0] *= omegaStar  # k
                temp_array[-1, 1] *= fStar**2   # P_phi
                temp_array[-1, 2] *= (fStar*omegaStar)**2  # P_phi'
                temp_array[-1, 4] *= omegaStar  # rho_k
                """
            else:  # if data is take at 'next time'
                if spectra is None:
                    spectra = np.array([temp_array])
                else:
                    spectra = np.append(spectra, np.array([temp_array]),
                                        axis=0)
                temp_array = None  # empty temp_array
    return spectra


def read_spec():
    """
    Read spectra data and interprete it
    """
    print("Reading spcetra...")
    data = np.genfromtxt("./average_spectra_times.txt").T  # read spectratimes
    # all in physical time
    spec_times = data
    spectra = parse_spec('./spectra_scalar_0.txt')
    return spec_times, spectra


def give_gradual_rgb(x):
    """
    Calculate color from a number
    x: number between 0 and 1
    return: RGB color (3-tuple with value between 0 and 1)
            range from red to yellow to blue according to x
    """
    if x <= 0.5 and x >= 0:  # from red to yellow
        return (1, 2*x, 0)
    elif x > 0.5 and x <= 1:  # from yellow to blue
        return (1-2*(x-0.5), 1-2*(x-0.5), 2*(x-0.5))
    else:
        print("Error in input of give_gradual_rgb()")
        return TypeError


def draw_spectra(k, spec, spec_t, range=None, k_range=None):
    """
    Draw spectra (one column from spectra.txt file)
    at different times in different colors.
    Arguments:
        spec is the big spec array from read_spectra.
        spec_t is the array from spectra_times.txt
        range is the range of time for plotting
    Return
        matplotlib Axes, ready to plot
    """
    if range is not None:  # if range for t is defined
        time_mask = (spec_t > range[0]) & (spec_t < range[1])
        spec = spec[time_mask]
        spec_t = spec_t[time_mask]

    if k_range is not None:  # if range for k is defined
        # make mask along 0th time, k-array is the same for different times
        mask = (k[0] >= k_range[0]) & (k[0] <= k_range[1])
        k = k[:, mask]
        spec = spec[:, mask]

    fig, ax = plt.subplots()
    for index, t in list(enumerate(spec_t)):
        color = give_gradual_rgb(float(index)/spec_t.shape[0])  # RGB tuple

        if index == 0 or index == spec_t.shape[0]-1:
            # to show the time range in the legend
            # put the first and last curve in the legend
            label = r"$t = %.2f / \omega_*$" % t
            ax.plot(k[index], spec[index], color=color, label=label)
        else:
            # else for everything else
            ax.plot(k[index], spec[index], color=color)

    return fig, ax


def plot_fluc_spec(spec_times, spectra, para, spec_range=None, k_range=None,
                   suffix=None, label=None):
    print("Plotting field fluctuation spectra...")
    k = spectra[:, :, 0] * para.omegaStar / para.m_eff_min
    fluc_spec = spectra[:, :,1]
    fig, ax = draw_spectra(k, fluc_spec, spec_times,
                           range=spec_range, k_range=k_range)
    ax.plot([], [], label=label)
    ax.set_xlabel(r"(com.) $k/|m_{\mathrm{eff, min}}|$")
    ax.set_ylabel(r"$P_{{\phi}}/\phi_0^2$")
    #  ax.set_yscale('log')
    ax.legend()

    plot_fn = "P_k_color"
    if k_range is not None or range is not None:
        plot_fn += "_zoomed_in"
    if suffix is not None:
        plot_fn += suffix
    plot_fn += ".pdf"
    plt.savefig(plot_fn, bbox_inches="tight")
    plt.close()


def find_sign_change(x):
    """
    Find sign changes of oscillation
    """
    sign = -1  # initially always negative
    sign_change_ind_array = np.array([0])

    for t_ind, t in list(enumerate(x)):
        #  t_ind_init = 0  # always starts with zero
        if x[t_ind] * sign > 0:  # x has the same sign as 'sign'
            #  t_ind_end = t_ind
            pass
        else:  # different signs
            sign_change_ind_array = np.append(sign_change_ind_array,
                                              int(t_ind))
            #  print("Sign change at t =", t, "H0")
            #  t_ind_init = t_ind
            sign *= -1  # flip sign
    return sign_change_ind_array


def plot_fluc_spec_by_vel(spec_times, spectra, average_scalar):
    """
    Plot fluctuation spectra according to sign changes in inflaton velocity
    """
    # check dimension of time array
    if spec_times.shape[0] != average_scalar[0].shape[0]:
        return ValueError

    fld_vld = average_scalar[2]
    sign_changes = find_sign_change(fld_vld)
    #  print(sign_changes)
    for i in range(0, sign_changes.shape[0]-1):
        times_ind = int(sign_changes[i])
        times_ind2 = int(sign_changes[i+1])
        t_range = [spec_times[times_ind], spec_times[times_ind2]]
        plot_fluc_spec(spec_times, spectra, spec_range=t_range,
                       label=r"$\dot{\phi} < 0$", suffix=str(i))
        #  print(times_ind, times_ind2, spec_times[times_ind], spec_times[times_ind2])
    return None


def find_max_P_k(k, spec, spec_t):
    """
    Find peak value of a spectrum
    """
    max_array = None
    for index, t in list(enumerate(spec_t)):
        # plot dots at maximum
        P_max = np.amax(spec[index])
        mask = (spec[index] == P_max)
        k_max = k[index, mask][0]  # after masking, still an array

        temp = np.array([[t, k_max, P_max]])
        if max_array is None:  # if not defined
            max_array = temp
        else:
            max_array = np.append(max_array, temp, axis=0)
    return max_array


def plot_max_P_k(spectra, spec_times, average_scalar, para):
    """
    Plot peak value of spectrum at different times
    """
    k = spectra[:, :, 0] / para.m_eff_min
    fluc_spec = spectra[:, :, 1]

    max_array = find_max_P_k(k, fluc_spec, spec_times).T
    t = max_array[0]
    k_max = max_array[1]
    P_max = max_array[2]

    fig, ax = plt.subplots()

    ax.plot(t, k_max, color="steelblue")
    ax.tick_params(axis="y", labelcolor="steelblue")
    ax.set_xlabel(r"$t\cdot \omega_*$")
    ax.set_ylabel(r"$k_{max}/|m_\mathrm{eff, min}|$", color="steelblue")
    ax.set_ylim((0, 2))
    ax.set_xlim((0.5, 2))

    ax2 = ax.twinx()
    ax2.plot(average_scalar[0], average_scalar[1]/para.phi0,
             linestyle="--", color="black")
    ax2.tick_params(axis="y", labelcolor="black")
    #  ax2.set_ylim(())
    ax2.set_ylabel(r"$\langle \phi \rangle / \phi_0$")

    plt.savefig("P_k_k_max.pdf", bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(t, P_max)
    ax.set_xlabel(r"$t\cdot \omega_*$")
    ax.set_ylabel(r"$P_{\phi, \mathrm{max}}/\phi_0^2$")
    plt.savefig("P_k_P_max.pdf", bbox_inches="tight")
    plt.close()


def plot_n_k(spec_times, spectra, para, spec_range=None, k_range=None):
    """
    Plot occupation number of produced particles; keep lattice unit of n_k
    """
    print("Plotting occupation number...")
    k = spectra[:, :, 0] / para.m_eff_min
    n_spec = spectra[:, :, 3]
    fig, ax = draw_spectra(k, n_spec, spec_times,
                           range=spec_range, k_range=k_range)
    ax.set_ylabel(r"$\tilde{n}_{\tilde{k}}$")
    ax.set_xlabel(r"(com.) $k/|m_{\mathrm{eff, min}}|$")
    #  ax2.set_yscale("log")
    ax.legend()
    if spec_range is None and k_range is None:
        plt.savefig("n_k_color.pdf", bbox_inches="tight")
    else:
        plt.savefig("n_k_color_zoomed_in.pdf", bbox_inches="tight")
    plt.close()


def plot_rho_k(spec_times, spectra, para, spec_range=None, k_range=None):
    """
    Plot Fourier transformed energy perturbations
    """
    print("Plotting energy density fluctuations...")
    k = spectra[:, :, 0] * para.omegaStar / para.m_eff_min
    rho_spec = spectra[:, :,4]
    fig, ax = draw_spectra(k, rho_spec, spec_times,
                           range=spec_range, k_range=k_range)
    ax.set_ylabel(r"$\delta \tilde{\rho}_{\tilde k}$")
    ax.set_xlabel(r"(com.) $k/|m_{\mathrm{eff, min}}|$")
    #  ax2.set_yscale("log")
    ax.legend()
    if spec_range is None and k_range is None:
        plt.savefig("rho_k_color.pdf", bbox_inches="tight")
    else:
        plt.savefig("rho_k_color_zoomed_in.pdf", bbox_inches="tight")
    plt.close()


def get_total_spec(spec_times, spectra, which_spec):
    """
    Sum over spec at one time instance
    """
    #  norm = (L * omegaStar)**3
    # keep everything in program untis
#      if which_spec == 3:  # n
    #      norm = omegaStar**3
    #  elif which_spec == 4:  # rho
    #      norm = omegaStar**4
    #  else:
    #      return ValueError
    spec_tot = np.zeros(spec_times.shape[0])
    for i in range(0, spectra.shape[0]):  # for every time
        temp_spec = spectra[i].T
        # multiply the number of points in the lattice
        n_k = temp_spec[which_spec] * temp_spec[-1]
        #  spec_tot[i] = np.sum(n_k) * norm
        spec_tot[i] = np.sum(n_k)
    return spec_tot


def plot_total_n(spec_times, spectra):
    """
    Plot total number of inflaton against time
    """
    print("Plotting total occupation number...")
    n_k_tot = get_total_spec(spec_times, spectra, 3)
    plt.plot(spec_times, n_k_tot, label="with CosmoLattice")
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.ylabel(r"$n_{\phi} / \omega_*^3 $")
    plt.yscale('log')
    plt.legend()
    plt.savefig("n_k_of_t.pdf", bbox_inches="tight")
    plt.close()


def plot_P_delta(spec_times, spectra, energy_array, para,
                 range=None, k_range=None, scale=None, suffix=None):
    """
    Plot power spectrum of density contrast
    """
    # TODO: something is wrong here...
    print("Plotting P_delta...")

    # appears in formula for \delta_k
    pre_fact = 1 / (2*np.pi**2)/para.L**3 / para.N**12
    # appears for different units of rho_k and average_rho
    pre_fact *= 1 / (para.omegaStar**2 * para.fStar**4)
    #  fig, ax = plt.subplots()

    if spec_times.shape[0] == energy_array[0].shape[0]:
        k_array = spectra[:, :, 0] * para.omegaStar
        P_delta_array = np.zeros(k_array.shape)
        for t_ind, t in enumerate(spec_times):
            # calculate P_delta
            aver_rho = energy_array[-1, t_ind]
            temp_spec = spectra[t_ind].T
            rho_k = temp_spec[-2]
            k = k_array[t_ind]
            P_delta = pre_fact * k**3 * (rho_k/aver_rho)**2

            P_delta_array[t_ind] = P_delta
    else:
        # frequencies for energies and spectra are different
        k_array = spectra[:, :, 0] * para.omegaStar
        P_delta_array = np.zeros(k_array.shape)
        t_aver_energy = energy_array[0]
        for t_ind, t in enumerate(spec_times):
            aver_index_tup = np.where(t_aver_energy == t)  # returns a tuple of array(s)
            if aver_index_tup[0].shape[0] != 1:
                print("Error in plot_P_delta()")
                return ValueError
            else:
                aver_index = aver_index_tup[0][0]
                #  print(t_aver_energy[aver_index], t)
                aver_rho = energy_array[-1, aver_index]
                temp_spec = spectra[t_ind].T
                rho_k = temp_spec[-2]
                k = k_array[t_ind]
                P_delta = pre_fact * k**3 * (rho_k/aver_rho)**2

                P_delta_array[t_ind] = P_delta
    # k_array already in planck units
    fig, ax = draw_spectra(k_array/para.m_eff_min, P_delta_array,
                           spec_times, range=range, k_range=k_range)
    ax.set_xlabel(r"(com.) $k/|m_{\mathrm{eff, min}}|$")
    ax.set_ylabel(r"$P_{\delta}$")
    if scale is not None:
        ax.set_yscale(scale)
    plt.legend()
    plot_fn = "P_delta"
    if k_range is not None or range is not None:
        plot_fn += "_zoomed_in"
    if suffix is not None:
        plot_fn += suffix
    plot_fn += ".pdf"
    plt.savefig(plot_fn, bbox_inches="tight")
    plt.close()
