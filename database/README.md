# **`Optab` opacity database**

This directory stores opacity databases used in Optab.

> ### Notes
> - Keep the directory structure as is because `Optab` programs use the relative paths.
> - As for molecular line lists, we recommend to start with HITRAN/HITEMP since Exomol database is gigantic.

---

### `1016620_Supplementary_Data/`
- This directory stores free-free Gaunt factor data by [van Hoof et al. (2014)](https://academic.oup.com/mnras/article/444/1/420/1016620).
- Extract their [supplementary data](https://academic.oup.com/mnras/article/444/1/420/1016620#supplementary-data) in this directory.

### `Karzas_Latter_1961.tsv`
- This file is "Table I" (bound-free Gaunt factors) in [Karzas and Latter (1961)](http://articles.adsabs.harvard.edu/pdf/1961ApJS....6..167K) (&copy; AAS. Reproduced with permission).

### `photo/`
- This directory contains [Verner's photoionization cross sections data](https://www.pa.uky.edu/~verner/photo.html).
- Copy all files in https://www.pa.uky.edu/~verner/dima/photo/ into this directory.

### `TOPbase/`
- This directory is a workspace for [TOPbase: photoionization cross sections](http://cdsweb.u-strasbg.fr/topbase/xsections.html).
1. Execute `get_topbase.sh` to retrieve the cross section data:
   ```
   % bash ../fetch/get_topbase.sh
   ```
2. Execute `convert_topbase_h5` to store the data in the specific HDF5 format for `Optab`:
   ```
   % ../src/convert_topbase_h5
   ```

### `NIST/`
- This directory is a workspace for NIST [Atomic Weights and Isotopic Compositions with Relative Atomic Masses](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses) and [Atomic Spectra Database Levels](https://physics.nist.gov/PhysRefData/ASD/levels_form.html).
1. Execte `get_nist_parallel.py` to retrieve the data:
   ```
   % python3 ../fetch/get_nist_parallel.sh
   ```
2. Execute `convert_nist_h5` to store the data in the specifice HDF5 format for `Optab`:
   ```
   % ../src/convert_nist_h5
   ```

### `HITRAN/`
- This directory is a workspace for [HITRAN](https://hitran.org/) molecular linelists.
1. Execulte `get_hitran_meta.sh` to retrieve [HITRAN Isotopologue Metadata](https://hitran.org/docs/iso-meta/):
   ```
   % bash ../fetch/get_hitran_meta.sh
   ```
> The following procedure is for H<sub>2</sub>O in HITRAN. Repeat it for other species.
2. Get linelists (`.par` files):
   1. Goto [Line-by-Line Search](https://hitran.org/lbl/).
   2. "Select Molecules" &rarr; check 1. H2O
   3. "Select Isotopologues" &rarr; check all isotopologues
   4. "Select Wavenumber / Wavelength Range" &rarr; leave blank for &nu;<sub>max</sub>
   5. "Select or Create Output Format" &rarr; do nothing
   6. "Start Data Search> Search Results" &rarr; download the "Output transitions data (160-character `.par` format)" as `original/01_HITRAN.par`. Here, `"01"` is the two digits [molecule ID](https://hitran.org/docs/molec-meta/) of H<sub>2</sub>O.
3. Break down the downloaded `.par` file to make separate `.par` files for different isotoplogue:
   ```bash
   $ ../src/preproc_hitran original/01_HITRAN.par
   ```
4. Convert each isotopologue `.par` file to an HDF5 file in the specific format for `Optab`:
   ```bash
   $ ../src/convert_lines_h5
   ```   
> The following procedure is for H<sub>2</sub>O in HITEMP. Repeat it for other species.
2. Get linelists (`.par` files):
   1. Goto https://hitran.org/hitemp/data/HITEMP-2010/H2O_line_list/.
   2. Download all zipped files and unzip them.
   3. Concatenate all unzipped files into a single file, `original/01_HITEMP2010.par`.

3. Break down the downloaded `.par` file to make separate `.par` files for different isotoplogue:
   ```bash
   $ ../src/preproc_hitran original/01_HITEMP2010.par
   ```
4. Convert each isotopologue `.par` file to an HDF5 file in the specific format for `Optab`:
   ```bash
   $ ../src/convert_lines_h5
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
1. Execute `get_kurucz_linelists.sh`, which retrieves two linelists, `gfall08oct17.dat` and `gfpred26apr18.dat`, from Kurucz database and puts them in the subdirectory `linelists/`. Edit `list_convert.txt` to list the relative paths of the linelists, and execute `convert_lines_h5` to generate HDF5 files for `Optab`:
   ```bash
   $ bash ../fetch/get_kurucz_linelists.sh
   $ ls linelists/*.dat > list_convert.txt
   $ ../src/convert_lines_h5
   ```
2. Execute `get_kurucz_gfgam.sh` to get the level data file `gf????.gam` for all species available, and execute `convert_gfgam_h5` to generate an HDF5 file for `Optab`:
   ```bash
   $ bash ../fetch/get_kurucz_gfgam.sh
   $ ../src/convert_gfgam_h5
   ```

### `h5/`
- This directory is the storage for the files in the HDF5 format made for `Optab`.
