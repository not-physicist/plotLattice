import numpy as np
import matplotlib.pyplot as plt

plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def read_average_scalar():
    """
    Read average field values and etc
    """
    data = np.genfromtxt("./average_scalar_0.txt").T
    # keep everything in program units
    #  data[1] *= fStar
    #  data[2] *= fStar * omegaStar
    #  data[3] *= fStar**2
    #  data[4] *= (fStar * omegaStar)**2
    #  data[5] *= fStar
    #  data[6] *= fStar * omegaStar
    return data


def plot_mean_fld(data, t_range=None):
    """
    Plot volume averaged field value
    """
    print("Plotting mean field value...")

    t = data[0]
    #  mean_fld = data[1] / phi0
    mean_fld = data[1]
    plt.plot(t, mean_fld)
    t0_array = np.array([t[0], t[-1]])
    #  plt.plot(t0_array, [1/3, 1/3], color="grey", linestyle="--", label="tachyonic")

    plt.fill_between(t0_array, 1/3, 1, color="lightgrey", label="tachyonic")

    plt.ylabel(r"$\langle \phi \rangle /\phi_0$")
    plt.legend()
    plt.xlabel(r"$t \cdot \omega_*$")
    fn = "mean_field"
    if t_range is None:
        plt.xlim((t[0], t[-1]))
    else:
        plt.xlim(t_range)
        fn += "_zoomed_in"
    plt.savefig(fn + ".pdf", bbox_inches="tight")
    plt.close()


def plot_rms(data):
    """
    Plot rms of inflaton field, i.e. sqrt(<phi^2> - <phi>^2)
    """
    print("Plotting field rms...")
    t = data[0]
    rms_fld = data[5]
    #  mask = (t > 0)  # only plot the interesting part
    #  plt.plot(t[mask], rms_fld[mask], label=r"$\sqrt{<\phi^2>}$")
    plt.plot(t, rms_fld)
    #  plt.legend()
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.ylabel(r"$\sqrt{\langle \phi^2 \rangle - \langle \phi \rangle^2}/\phi_0$")
    plt.savefig("rms.pdf", bbox_inches="tight")
    plt.close()

    rms_fld_vlc = data[6]
    plt.plot(t, rms_fld_vlc)
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.ylabel(r"$\sqrt{\langle \phi'^2 \rangle - \langle \phi' \rangle^2}/\phi_0/omega_*$")
    plt.savefig("rms_vlc.pdf", bbox_inches="tight")
    plt.close()


def plot_phi2(data):
    """
    Plot two types of phi^2
    """
    print("Plotting phi^2...")
    t = data[0]
    mean_fld = data[1]
    mean_fld2 = data[3]
    plt.plot(t, mean_fld2, label=r"$<\phi^2>$")
    plt.plot(t, mean_fld**2, label=r"$<\phi>^2$", color="black", linestyle=":")
    plt.legend()
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.ylabel(r"$\phi^2/\phi_0^2$")
    plt.savefig("phi2.pdf", bbox_inches="tight")
    plt.close()


def plot_fld_vel(data):
    """
    Plot field velocity, i.e. phi'
    """
    print("Plotting field velocities...")
    t = data[0]
    mean_vel = data[2]
    plt.plot(t, mean_vel)
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.ylabel(r"$\langle \dot{\phi} \rangle /(\phi_0 \omega_*)$")
    plt.savefig("phi_dot.pdf", bbox_inches="tight")
    plt.close()


def read_scale():
    """
    Read data for scale factors
    """
    print("Reading scale factors...")
    data = np.genfromtxt("./average_scale_factor.txt").T
    # keep everything in program units
    #  data[0] *= 1/omegaStar * H0
    #  data[2] *= omegaStar
    #  data[3] *= omegaStar
    return data


def plot_scale(data):
    """
    Plot scale factor against time
    """
    print("Plotting scale factors...")
    t = data[0]
    scale_fact = data[1]
    plt.plot(t, scale_fact)
    plt.yscale("log")
    #  plt.legend()
    plt.ylabel(r"$a$")
    plt.xlabel(r"$t \cdot \omega_*$")
    plt.savefig("scale_factor.pdf", bbox_inches="tight")
    plt.close()


def read_energies():
    data = np.genfromtxt("./average_energies.txt").T
    # keep everything in program units
    #  data[0] *= 1/omegaStar * H0
    #  data[1:] *= (fStar*omegaStar)**2
    return data


def get_w(data):
    """
    Calculate equation of state parameter omega from energy array;
    a single scalar field, three term in potential
    """
    t = data[0]
    K = data[1]
    G = data[2]
    V_tot = data[3] + data[4] + data[5]
    rho_tot = data[6]
    p_tot = K - G/3 - V_tot
    # no need to convert unit, cancelled anyway
    w = p_tot/rho_tot

    return t, w


def plot_w(energy_array):
    """
    Plot equation of state parameter omega
    """
    t, w = get_w(energy_array)
    fig, ax = plt.subplots()
    ax.plot(t, w)
    ax.set_ylabel(r"$\omega$")
    ax.set_xlabel(r"$t \cdot \omega_*$")
    plt.savefig("EOS-w.pdf", bbox_inches="tight")
    plt.close()


def draw_energies(energy_array, scale_fact):
    """
    Draw different comoving energies: V, G, K
    """
    if energy_array[0].shape[0] != scale_fact[0].shape[0]:
        print("ERROR in draw_energies()!")
        return ValueError
    t = energy_array[0]
    rho = energy_array[-1]
    G = energy_array[2]
    K = energy_array[1]
    V = (energy_array[3]+energy_array[4]+energy_array[5])
    # convert to per comoving volume
    a = scale_fact[1]
    com_G = G * a**3
    com_rho = rho * a**3
    com_K = K * a**3
    com_V = V * a**3

    fig, ax = plt.subplots()
    ax.plot(t, com_rho, label=r"$\rho_{tot}$")
    ax.plot(t, com_G, label=r"$\rho_G$")
    ax.plot(t, com_K, label=r"$\rho_K$")
    ax.plot(t, com_V, label=r"$\rho_V$")

    return fig, ax


def plot_energies(energy_array, scale_fact):
    """
    Plot different comoving energies: V, G, K
    """
    fig, ax = draw_energies(energy_array, scale_fact)
    ax.set_xlabel(r"$t \cdot \omega_*$")
    ax.set_ylabel(r"$\rho_i a^3 / (\phi_0^2 \omega_*^2)$")
    ax.legend()
    plt.savefig("comoving_energies.pdf", bbox_inches="tight")
    plt.close()

