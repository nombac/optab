import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
import argparse

def plot_n_species(fname, syms):
    """
    Plot number density for each species on the T-P plane.
    """
    with h5py.File(fname, 'r') as h5file:
        n_layer = h5file['n_layer'][()]
        n_species = h5file['n_species'][()]
        ndens = h5file['ndens'][()]
        temp = h5file['temp'][()]
        pres = h5file['pres'][()]
        species = h5file['species'][()]

    # Decode species names if stored as bytes
    species_names = [s.decode() if isinstance(s, bytes) else s for s in species]

    # Create output directory
    outdir = fname.replace('.h5', '')
    os.makedirs(outdir, exist_ok=True)

    xrange = [1e2, 1e8]
    yrange = [1e-5, 1e10]
    vrange = [0, 25]

    cmap = plt.cm.jet

    for i in range(n_species.item()):
        title = species_names[i].strip()
        if len(title) == 0:
            continue

        nd = ndens[:, i]

        # Check for NaN/Inf
        if not np.all(np.isfinite(nd)):
            print(f"Warning: non-finite values in {title}, layer indices: {np.where(~np.isfinite(nd))[0]}")
            continue

        if np.max(nd) == 0.0:
            continue

        v = np.log10(nd, where=(nd > 0), out=np.full_like(nd, dtype=float, fill_value=vrange[0]))

        norm = mcolors.Normalize(vmin=vrange[0], vmax=vrange[1])

        fig, ax = plt.subplots()
        sc = ax.scatter(temp, pres, c=v, s=syms, marker='s', cmap=cmap, norm=norm)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_xlabel('T [K]')
        ax.set_ylabel('P [Ba]')
        ax.set_title(title)

        cbar = plt.colorbar(sc)
        cbar.set_label(r'log n [cm$^{-3}$]', rotation=270, labelpad=15)

        ax.text(0.8, 0.95, outdir, transform=ax.transAxes,
                fontsize=6, color='gray', fontfamily='monospace')

        png_name = os.path.join(outdir, title + '.png')
        plt.savefig(png_name, format='png', transparent=True, dpi=600)
        plt.close(fig)
        print(png_name)

def main():
    parser = argparse.ArgumentParser(description='Plot number density for each species from an HDF5 file.')
    parser.add_argument('filename', type=str, help='Name of the HDF5 file.')
    parser.add_argument('--syms', type=int, default=5, help='Symbol size for the scatter plot (default: 5).')
    args = parser.parse_args()

    plot_n_species(args.filename, args.syms)

if __name__ == "__main__":
    main()
