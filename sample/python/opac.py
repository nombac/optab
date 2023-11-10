import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse

# Constants
K_BOL = 1.3806503e-16

def read_mono_hdf5(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            return {name: f[name][0] for name in f.keys()}
    except Exception as e:
        raise RuntimeError(f"Failed to read {file_path}: {e}")

def opac(dir_path, mean, syms):
    files = [os.path.join(dir_path, 'output', f) for f in os.listdir(os.path.join(dir_path, 'output')) if f.startswith("mono_") and f.endswith(".h5")]
    nmax = len(files)
    
    if nmax == 0:
        raise ValueError("No HDF5 files found.")

    print(f"Reading {nmax} layers ... ")

    datasets = [read_mono_hdf5(file) for file in files]
    tmp2 = np.array([data['temp2'] for data in datasets])
    tmp_tot = np.array([data['temp'] for data in datasets])
    rho_tot = np.array([data['rho'] for data in datasets])
    pre_tot = np.array([data['nden'] * K_BOL * data['temp'] for data in datasets])
    ros_tot = np.log10(np.array([data['ros'] for data in datasets]) / rho_tot)
    pla_tot = np.log10(np.array([data['plac'] + data['plal'] for data in datasets]) / rho_tot)
    pla2_tot = np.log10(np.array([data['plac2'] + data['plal2'] for data in datasets]) / rho_tot)

    # Color scaling for visualization
    vmax, vmin = 7.0, -6.0

    # Title and label mapping
    titles = {
        'ross': ('Rosseland-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
        'pla': ('Planck-mean opacity', 'log $\kappa$ [cm$^2$/g]'),
        'pla2': (f'Planck-mean opacity at T_rad={int(tmp2[0])}K', 'log $\kappa$ [cm$^2$/g]')
    }

    if mean not in titles:
        raise ValueError(f"Unknown mean type: {mean}")
    title, btitle = titles[mean]

    # Plotting
    plt.rcParams.update({'font.size': 16})  # Adjust the number to your preference
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(tmp_tot, rho_tot, c=
                          ros_tot if mean == 'ross' else
                          pla_tot if mean == 'pla' else
                          pla2_tot,
                          s=syms, cmap='jet', norm=plt.Normalize(vmin=vmin, vmax=vmax), marker='s')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('T [K]')
    plt.ylabel(r'$\rho$ [g/cm$^3$]')
    plt.title(title)
    plt.colorbar(scatter, label=btitle)
    plt.grid(True)

    # Save the plot in different formats
    file_basename = os.path.join(dir_path, 'output', mean)
    for ext in ['pdf', 'png']:
        plt.savefig(f'{file_basename}.{ext}', format=ext, transparent=True)
        print(f'Output file: {file_basename}.{ext}')

    plt.show()

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process some integers.')
    # Add the arguments
    parser.add_argument('dir_path', type=str, help='The path to the directory containing HDF5 files')
    parser.add_argument('mean', type=str, choices=['ross', 'pla', 'pla2'], help='Type of mean opacity to plot')
    parser.add_argument('syms', type=int, help='Marker size for scatter plot')

    # Execute the parse_args() method
    args = parser.parse_args()

    opac(args.dir_path, args.mean, args.syms)

if __name__ == '__main__':
    main()
