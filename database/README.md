# **`Optab` opacity database**

This directory stores opacity databases used in Optab.

> ### Notes
> - Keep the directory structure as is because `Optab` programs use the relative paths.
> - As for molecular line lists, we recommend to start with HITRAN/HITEMP since Exomol database is gigantic.

---
### `h5/`
- This directory is the storage for the files in the HDF5 format made for `Optab`.
- Make this directory first.
   ```
   mkdir h5/
   ```

### `1016620_Supplementary_Data/`
- This directory stores free-free Gaunt factor data by [van Hoof et al. (2014)](https://academic.oup.com/mnras/article/444/1/420/1016620).
- Download their [supplementary data](https://academic.oup.com/mnras/article/444/1/420/1016620#supplementary-data) and extract it in this directory.
   ```
   unzip ~/Downloads/1016620_Supplementary_Data.zip -d 1016620_Supplementary_Data
   ```

### `Karzas_Latter_1961.tsv`
- This file is "Table I" (bound-free Gaunt factors) in [Karzas and Latter (1961)](http://articles.adsabs.harvard.edu/pdf/1961ApJS....6..167K) (&copy; AAS. Reproduced with permission).

### `photo/`
- This directory stores [Verner's photoionization cross sections data](https://www.pa.uky.edu/~verner/photo.html).
- Copy all files in https://www.pa.uky.edu/~verner/dima/photo/ into this directory.
   ```
   wget -r -np -nH --cut-dirs=3 -P photo -R "index.html*" https://www.pa.uky.edu/~verner/dima/photo/
   ```

### `TOPbase/`
- This directory is a workspace for [TOPbase: photoionization cross sections](http://cdsweb.u-strasbg.fr/topbase/xsections.html).
1. Execute `get_topbase.py` to retrieve the cross section data files and convert them to a specific HDF5 format for `Optab`:
   ```
   cd TOPbase/
   ```
   ```
   python3 ../fetch/get_topbase.py
   ```

### `NIST/`
- This directory is a workspace for NIST [Atomic Weights and Isotopic Compositions with Relative Atomic Masses](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses) and [Atomic Spectra Database Levels](https://physics.nist.gov/PhysRefData/ASD/levels_form.html).
1. Execute `get_nist_parallel.py` to retrieve the level/atomic data and convert them to a specific HDF5 format for `Optab`:
   ```
   cd NIST/
   ```
   ```
   python3 ../fetch/get_nist_parallel.py
   ```
   or try `../fetch/get_nist.py` (slower) if you encounter a network issue.

### `HITRAN/`
- This directory is a workspace for [HITRAN](https://hitran.org/) molecular linelists.
1. Execute `get_hitran.py` to retrieve [HITRAN Isotopologue Metadata](https://hitran.org/docs/iso-meta/) and the partition function files:
   ```
   cd HITRAN/
   ```
   ```
   python3 ../fetch/get_hitran_Qs.py
   ```
1. Get linelists (`.par` files):
   1. `HITRAN` lines (e.g. H2O; repeat this procedure for other species.)
      1. Goto [Line-by-Line Search](https://hitran.org/lbl/).
      1. "Select Molecules" &rarr; check 1. H2O
      1. "Select Isotopologues" &rarr; check all isotopologues
      1. "Select Wavenumber / Wavelength Range" &rarr; leave blank for &nu;<sub>max</sub>
      1. "Select or Create Output Format" &rarr; .par (160 chars)
      1. "Start Data Search> Search Results" &rarr; download the "Output transitions data (160-character `.par` format)" as `original/01_HITRAN.par`. Here, `"01"` is the two digits [molecule ID](https://hitran.org/docs/molec-meta/) of H<sub>2</sub>O.
   1. `HITEMP` lines
      1. Go to [HITEMP](https://hitran.org/hitemp/) and download bzip2ed `.par` files to `original/` and bunzip2 them. Note that tne data for H<sub>2</sub>O and CO<sub>2</sub> is divided into multiple files, respectively. In these cases, edit `../fetch/get_hitemp_multi.sh` appropriately and run it to get a single .par file:
         ```bash
         bash ../fetch/get_hitemp_multi.sh
         ```
1. Break down the downloaded `.par` files in `original/` to make separate `.par` files for different isotoplogue, and convert them to HDF5 files:
   ```bash
   bash ./preproc_and_convert_HITRAN.sh
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

### `Kurucz/`
- This directory is a workspace for [Kurucz atomic linelinsts](http://kurucz.harvard.edu/linelists.html)
1. Execute `get_kurucz_linelists.sh` to retrieve two linelists, `gfall08oct17.dat` and `gfpred26apr18.dat`, from Kurucz database and convert them to HDF5 files for `Optab`:
   ```bash
   bash ../fetch/get_kurucz_linelists.sh
   ```
2. Execute `get_kurucz_gfgam.sh` to get the level data files `gf????.gam` for all species available (ignore `Not Found` errors) and convert them to an HDF5 file for `Optab`:
   ```bash
   python3 ../fetch/get_kurucz_gfgam.py
   ```
