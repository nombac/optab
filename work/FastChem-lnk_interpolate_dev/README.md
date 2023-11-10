# Calculating a chemical abundance table with FastChem
Here, by adopting the branch `lnk_interpolate_dev`, you run `FastChem` with tabulated lnK files, not with fitting functions of lnK. A program named `prep_FastChem` makes lnK files for molecules from the original fitting functions while those for atomic ions from the energy level data used in `Optab`.

**WARNING:** For some atomic ions at high temperatures, the values of lnK may be inaccurate (underestimated) because they are computed from the NIST energy level data, which is not complete for higher levels.

1. Build `FastChem`
    1. `% cd <some_working_dir>`
    1. `% git clone https://github.com/exoclime/FastChem.git`
    1. `% cd FastChem/`
    1. `% git checkout lnk_interpolate_dev`
    1. Comment out the following line 37 in `fastchem_src/check.cpp`:\
    `if (this->number_density < min_limit) this->number_density = min_limit;`
    1. `% mkdir build`
    1. `% cd build/`
    1. `% cmake .. && make`
    1. `% cd ../`
1. Generate an initial configuration file (and its associated files) for `FastChem`
    1. `% cd input/`
    1. `% rsync -avLu $OPTAB/work/FastChem-lnk_interpolate_dev/input/ .`
    1. Change appropriately a preprocessor macro, `OPTAB_DATABASE_DIR`, defined at the top of `prep_FastChem.F90`.
    1. `% make prep_FastChem`
    1. Review `prep_FastChem.dat`
        1. label
        1. temperature grid
        1. pressure grid
    1. `% ./prep_FastChem`
    1. `% cd ../`
1. Run `FastChem`
    1. `% mkdir output/`
    1. `% ./fastchem input/init.config_<label>`