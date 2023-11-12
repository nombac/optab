# **`Optab` opacity database**

This directory houses the opacity databases utilized by Optab.

> ### Notes
> - Retain the current directory structure, as Optab relies on relative path references.
> - We recommend starting with the HITRAN/HITEMP database for molecular line lists, as the Exomol database's extensive size may be overwhelming for initial exploration.
> - Please be aware that the scripts provided for downloading the databases may become obsolete if there are changes in the database formats. If you encounter such issues, kindly inform the author (shirose@jamstec.go.jp).
---
### `h5/`
- This directory serves as the repository for HDF5-formatted files created for `optab``.

### `1016620_Supplementary_Data/`
- This directory contains the free-free Gaunt factor data authored by [van Hoof et al. (2014)](https://academic.oup.com/mnras/article/444/1/420/1016620).
1. Download their [supplementary data](https://academic.oup.com/mnras/article/444/1/420/1016620#supplementary-data) and extract it in this directory.
   ```bash
   unzip ~/Downloads/1016620_Supplementary_Data.zip -d 1016620_Supplementary_Data
   ```

### `Karzas_Latter_1961.tsv`
- This file is "Table I" (bound-free Gaunt factors) in [Karzas and Latter (1961)](http://articles.adsabs.harvard.edu/pdf/1961ApJS....6..167K) (&copy; AAS. Reproduced with permission).

### `photo/`
- This directory stores [Verner's photoionization cross sections data](https://www.pa.uky.edu/~verner/photo.html).
1. Copy all files in https://www.pa.uky.edu/~verner/dima/photo/ into this directory.
   ```bash
   wget -r -np -nH --cut-dirs=3 -P photo -R "index.html*" https://www.pa.uky.edu/~verner/dima/photo/
   ```

### `TOPbase/`
- This directory serves as a workspace for [TOPbase: photoionization cross sections](http://cdsweb.u-strasbg.fr/topbase/xsections.html).
   ```bash
   cd TOPbase/
   ```
1. Execute `get_topbase.py` to retrieve the cross section data files and convert them to a specific HDF5 format for `Optab`:
   ```bash
   python3 get_topbase.py
   ```

### `NIST/`
- This directory functions as a workspace for the NIST databases on [Atomic Weights and Isotopic Compositions with Relative Atomic Masses](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses) and [Atomic Spectra Database Levels](https://physics.nist.gov/PhysRefData/ASD/levels_form.html).
   ```bash
   cd NIST/
   ```
1. Execute `get_nist_parallel.py` to retrieve the level/atomic data and convert them to a specific HDF5 format for `Optab` (**REQUIREMENT: [`lynx`](https://lynx.invisible-island.net/)**):
   ```bash
   python3 get_nist_parallel.py
   ```
   or try `get_nist.py` (slower) if you encounter a network issue.

### `HITRAN/`
- This directory is a workspace for [HITRAN](https://hitran.org/) molecular linelists.
   ```bash
   cd HITRAN/
   ```
1. [HITRAN Isotopologue Metadata](https://hitran.org/docs/iso-meta/) and the partition function files (**REQUIREMENT: [`w3m`](https://w3m.sourceforge.net/)**):
   ```bash
   python3 get_hitran_meta.py
   ```
1. Linelists (`.par` files):
   1. [`LBL`](https://hitran.org/lbl/) (**REQUIREMENT: Goggle Chrome**):
      
      ```bash
      python3 get_hitran_lines.py "/Users/shirose/Library/Application\ Support/Google/Chrome/Default"
      ```
      The required argument is your Chrome user profile directory. To find this, visit chrome://version in Chrome..

      >NOTE:
      >- If you are redirected to the registration page, register or log in. Then, exit Chrome and restart the process.
      >- Alternatively, you can download the linelists manually; refer to the instructions provided within the code.
   1. [`HITEMP`](https://hitran.org/hitemp/)
      ```bash
      bash get_hitemp.sh
      bash get_hitemp_multi.sh
      ```
1. Create individual .par files for each isotopologue and convert them to HDF5 files ready for optab:
   ```bash
   bash preproc_and_convert_HITRAN.sh
   ```

### `Kurucz/`
- This directory serves as a workspace for [Kurucz atomic linelinsts](http://kurucz.harvard.edu/linelists.html)
   ```bash
   cd Kurucz/
   ```
1. Execute `get_kurucz_linelists.sh` to retrieve two linelists, `gfall08oct17.dat` and `gfpred26apr18.dat`, from Kurucz database and convert them to HDF5 files for `Optab`:
   ```bash
   bash get_kurucz_linelists.sh
   ```
2. Execute `get_kurucz_gfgam.sh` to get the level data files `gf????.gam` for all species available (ignore `Not Found` errors) and convert them to an HDF5 file for `Optab`:
   ```bash
   python3 get_kurucz_gfgam.py
   ```

### `Exomol/`
- This directory is a workspace for [Exomol](https://www.exomol.com/) molecular linelists. 
> The following procedure is for `1H2-16O__POKAZATEL`. Repeat it for other species.
1. Get a set of data files (`.def`, `.pf`, `.states`, `.trans` files) from https://www.exomol.com/data/molecules/H2O/1H2-16O/POKAZATEL/ and place them in the subdirectory `1H2-16O__POKAZATEL/`.
2. Get `.broad` files from https://www.exomol.com/data/molecules/H2O/, and place them in this directory. Note that the `.broad` files are provided only for limited species; if not provided, ignore this step.
3. Here is the list of the data files:
   ```
   1H2-16O__H2.broad
   1H2-16O__He.broad
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL.def
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL.pf
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL.states
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00000-00100.trans
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00100-00200.trans
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00200-00300.trans
   ...
   ```
4. Edit `list_convert.txt` to list the relative paths of the `.trans` files, and convert them into a n HDF5 file in the specific format for `Optab`:
   ```bash
   $ cat list_convert.txt
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00000-00100.trans
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00100-00200.trans
   1H2-16O__POKAZATEL/1H2-16O__POKAZATEL__00200-00300.trans
   ...
   $ ../src/convert_lines_h5
   ```
