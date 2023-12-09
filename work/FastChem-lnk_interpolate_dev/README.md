# Calculating a Chemical Abundance Table Using FastChem

In this section, we outline the process of running FastChem with tabulated lnK files, as opposed to using fitting functions for lnK. This is achieved by employing the `lnk_interpolate_dev` branch. The utility `prep_FastChem` is employed to generate lnK files for molecules. These files are derived from the original fitting functions. In contrast, lnK files for atomic ions are produced from the energy level data utilized in Optab.

**Warning**: It is important to note that for certain atomic ions at elevated temperatures, the lnK values may be subject to inaccuracies, specifically underestimations. This is because these values are calculated using NIST energy level data, which does not comprehensively cover higher levels.

---

1. **Build `FastChem`:**
   ```
   cd <some_working_dir>
   git clone https://github.com/exoclime/FastChem.git
   cd FastChem/
   git checkout lnk_interpolate_dev
   ```
   - Comment out the following line 37 in `fastchem_src/check.cpp`:\
   `if (this->number_density < min_limit) this->number_density = min_limit;`   
   ```
   mkdir build
   cd build/
   cmake .. && make
   cd ../
   ```
1. **Generate necessary files for `FastChem` including lnK files in `input/`:**
   ```
   cd $OPTAB/work/FastChem-lnk_interpolate_dev/input
   ```
   
   - Update the `OPTAB_DATABASE_DIR` preprocessor macro in `prep_FastChem.F90` to the correct path.
   - Update the `FASTCHEM_INPUT_DIR` preprocessor macro in `prep_FastChem.F90` to the correct path.

   ```
   make prep_FastChem
   ```
   
   - Verify the settings in` prep_FastChem.dat`, including the label, temperature grid, and pressure grid.
   
   ```
   ./prep_FastChem
   ```
1. **Run `FastChem:**
   ```
   cd $FASTCHEM
   ./fastchem input/config.input_<label>
   ```
   - Find the created chemical abundance table at `output/<label>.dat`. Here, `<label>` is `table` in the default setting.
