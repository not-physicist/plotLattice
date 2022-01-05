import h5py
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
#  import sys


def traverse_datasets(hdf_file):
    '''
    Go through a hdf file
    from: https://stackoverflow.com/questions/51548551/reading-nested-h5-group-into-numpy-array
    '''
    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = f'{prefix}/{key}'
            if isinstance(item, h5py.Dataset):  # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group):  # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    for path, _ in h5py_dataset_iterator(hdf_file):
        yield path


def list_cont(f, fn):
    '''
    List content of hdf5 file
    '''
    print(f"Content of {fn}:")
    for path in traverse_datasets(f):
        print(f"Path: {path}\tShape: {f[path].shape}\tData type: {f[path].dtype}")


def find_last(f):
    '''
    Find last entry of snapshots
    '''
    max_t = 0
    max_path = ''
    for path in traverse_datasets(f):
        # path name always is like /E_i/200.
        # choose part after / and delete the dot
        if path.endswith("."):
            s = path.split('/')[-1] + "0"
        else:
            s = path.split('/')[-1]
        t = int(float(s))
        if t > max_t:
            max_t = t
            max_path = path
    print(f"Last entry is {max_path}")
    return t, max_path


def get_CL_mean(t, which_E='tot'):
    '''
    Read CL's mean from average_energies.txt at time t
    '''
    index = -1
    if which_E != 'tot':
        # TODO: maybe implement this for other energy component;
        # right now not important
        print("Function in get_CL_mean() not implemented yet!")
        return ValueError

    CL_mean_fn = "./average_energies.txt"
    data = np.genfromtxt(CL_mean_fn)
    for line in data:
        #  print(np.abs(line[0]-t))
        if np.abs(line[0] - t) < 1e-3:
            #  print(line)
            #  print("found it")
            E = line[index]

    try:
        return E
    except Exception:
        print("Error in get_CL_mean()!")
        return ValueError


def check_cons(f, path, t, E_CL):
    '''
    Checking consistency:
        compute average once again and compare to CL's average
    '''
    if E_CL is not ValueError:
        dset = f[path]
        #  print(dset[:, :, :])
        flat_dset = dset[:, :, :].flatten()

        # consistency check
        cons = np.average(flat_dset)/E_CL
        if np.abs(cons - 1) < 1e-5:
            print(f"Consistency test passed for {path}")
            return True
        else:
            print(f"Consistency test failed for {path}")
            return False
    else:
        return ValueError


def find_num_over_thrhld(f, path, threshold):
    '''
    Find number of points in lattice which have more than threhold value
    '''
    dset = f[path]
    #  print(dset[:, :, :])
    flat_dset = dset[:, :, :].flatten()
    orig_num = flat_dset.shape[0]

    mask = (flat_dset > threshold)
    num = flat_dset[mask].shape[0]
    perc = 100*num/orig_num
    print(f"{num} of {orig_num} are above threshold. That is {perc}%")


def plot_histo(f, path, E_CL):
    '''
    Plot histogram of snapshots of energy components
    '''
    dset = f[path]
    #  print(dset[:, :, :])
    flat_dset = dset[:, :, :].flatten()

    mask = (flat_dset > 0)
    flat_dset = flat_dset[mask]
    ratio = flat_dset / E_CL

    bins = [i for i in range(0, 60)]
    bins.append(100)  # overflow
    #  hist, hist_edges = np.histogram(flat_dset, density=True, bins=bins)
    plt.hist(ratio, bins=bins)

    caption = f'Total number of points: \n{flat_dset.shape[0]: 1.2e}'
    plt.annotate(caption, xy=(0.65, 0.90), xycoords='axes fraction')
    plt.yscale("log")
    plt.ylabel("#")
    plt.xlabel(r"$\rho/\langle \rho \rangle$")
    plt.savefig("over_threshold_histo.pdf", bbox_inches="tight")
    plt.close()


def plot_over_threshold(f, path, threshold):
    '''
    Plot 3D scattering of points over threshold
    '''
    dset = f[path]

    # 3d scattering
    # creating figures
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    N = dset.shape[0]  # = # of grid points in lattice
    for i in range(0, N):
        for j in range(0, N):
            for k in range(0, N):
                E = dset[i, j, k]
                if E > threshold:
                    ax.scatter(i, j, k, color="red")
                elif E > threshold*0.8:
                    ax.scatter(i, j, k, color="orange")

    plt.show()
    #  plt.savefig("snapshots.png", bbox_inches="tight")
    #  plt.close()


def plot_2d_slice(f, path, CL_mean, z):
    '''
    z for which slice in z-direction
    '''
    dset = f[path]

    N = dset.shape[0]
    heatmap = dset[:, :, z] / CL_mean
    extent = [0, N-1, 0, N-1]

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    plt.clf()
    im = plt.imshow(heatmap, extent=extent)

    # colorbar
    cax = plt.axes([0.75, 0.1, 0.04, 0.8])  # [left, bottom, width, height]
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(r"$\rho/\langle \rho \rangle$")
    plt.xlabel("x / (L/N)")
    plt.xlabel("y / (L/N)")
    plt.savefig(f"./snap_2d_z={z}.pdf", bbox_inches="tight")

'''
#  plot_snap_over_threshold(200, 10)
fn = "./total_energy_snapshot_scalar.h5"
with h5py.File(fn, 'r') as f:
    #  list_cont(f, fn)
    t_last, path_last = find_last(f)
    CL_mean = get_CL_mean(t_last)
    #  check_cons(f, path_last, t_last, CL_mean)
    #  find_num_over_thrhld(f, path_last, 2*CL_mean)
    #  plot_histo(f, path_last, CL_mean)
    #  plot_over_threshold(f, path_last, 50*CL_mean)
    plot_2d_slice(f, path_last, CL_mean, 100)
'''
