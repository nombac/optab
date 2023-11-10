import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import subprocess

# Constants
K_BOL = 1.3806503e-16
C2 = 1.4387766844642551

def read_mono_hdf5(file_path):
    with h5py.File(file_path, 'r') as f:
        return {name: f[name][:] for name in f.keys()}

def generate_plot(dir_path, layer):
    files = sorted([f for f in os.listdir(os.path.join(dir_path, 'output')) if f.startswith("mono_") and f.endswith(".h5")])
    nmax = len(files)
    
    if layer > nmax - 1:
        raise ValueError(f'ERROR: layer {layer} > layer max {nmax-1}')

    print(f'Reading {nmax} layers ... ')
    data = read_mono_hdf5(os.path.join(dir_path, 'output', f'mono_{layer:05d}.h5'))
    
    tmp = data['temp']
    grd = data['grd']
    abs_ = data['abs']
    sca = data['sca']
    cnt = data['cnt']
    np_ = data['nden']
    pre = np_[0] * K_BOL * tmp

    grd0 = np.logspace(1, 10, 101)
    planck = grd0**3 / (np.exp(C2 * grd0 / tmp[0]) - 1)
    planck /= np.max(planck)
    
    ps_name = os.path.join(dir_path, 'output', f'mono_{layer:05d}')

    plt.rcParams.update({'font.size': 16})  # Adjust the number to your preference
    plt.figure(figsize=(10, 8))
    plt.loglog(grd, abs_, color='red', label='abs (cnt. + line)')
    plt.loglog(grd, cnt - sca, color='gray', label='abs (cnt.)')
    plt.loglog(grd, sca, color='blue', label='sca')
    plt.loglog(grd0, planck, color='green', linestyle=':', label='Planck')
    plt.xlabel(r'$\nu$ [cm$^{-1}$]')
    plt.ylabel(r'$\alpha$ [cm$^{-1}$]')
    plt.title('Monochromatic Opacity')
    plt.legend(loc='best')
    plt.ylim([1e-35,1e10])
    plt.grid(True)
    # Add text for tmp and pre
    plt.text(0.05, 0.15, f'T = {tmp.item():.2f} K', transform=plt.gca().transAxes, fontsize=12)
    plt.text(0.05, 0.10, f'P = {pre.item():.2e} Ba', transform=plt.gca().transAxes, fontsize=12)
    
    plt.savefig(f'{ps_name}.pdf', format='pdf', transparent=True)
    plt.savefig(f'{ps_name}.png', format='png', transparent=True)

    plt.show()
    print(f'Output files: {ps_name}.pdf, {ps_name}.png')

def main():
    parser = argparse.ArgumentParser(description='Monochromatic Opacity Plotter')
    parser.add_argument('dir_path', type=str, help='The path to the directory containing HDF5 files')
    parser.add_argument('layer', type=int, help='Layer number to process and plot')
    args = parser.parse_args()

    generate_plot(args.dir_path, args.layer)

if __name__ == '__main__':
    main()
