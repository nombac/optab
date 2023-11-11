# Preparing chemical abundance tables for `optab`

## Overview
A chemical abundance table consists of sets of *P* (pressure), *T* (temperature), and *N<sub>i</sub>* (number density of species *i*) like below, which should be converted to a specific HDF5 format to be read by `optab`.

   | P [bar] | T[K] | N<sub>H</sub> [cm<sup>-3</sup>] | N<sub>e</sub> [cm<sup>-3</sup>] | ...  |
   | ---------------- | ---------------- | ---------------- | ---------------- | ---- |
   | 1.0000000000e-12 | 2.5200000000e+03 | 2.6461950252e+06 | 2.9824566824e+02 | ...  |
   | 3.1622776602e-11 | 2.5200000000e+03 | 8.3691983426e+07 | 8.4813733226e+03 | ...  |
   | ...              | ...              | ...              | ...              | ...  |

## HDF5 format for `optab`

The HDF5 file of the chemical equilibrium abundances should have the following structure:
* `eos.h5`
  * `n_layer` (no. of dim = 1, dim. size = 1, data type = 32bit integer): number of layers
  * `n_species` (no. of dim = 1, dim. size = 1, data type = 32bit integer): dimension of the species array = 11000
  * `temp` (no. of dim = 1, dim. size = `n_layer`, data type = 64-bit floating-point): temperature
  * `pres` (no. of dim = 1, dim. size = `n_layer`, data type = 64-bit floating-point): pressure
  * `rho` (no. of dim = 1, dim. size = `n_layer`, data type = 64-bit floating-point): mass density (optional)
  * `ndens` (no. of dim = 2, dim. size = `n_layer` x `n_species`, data type = 64-bit floating-point): number density of species

The above file can be created by the following piece of Fortran90 code with HDF5 library (see for example, [convert_FastChem.F90](https://github.com/nombac/optab/blob/main/eos/src/convert_FastChem.F90). Here, `jmax` is the number of layers and `NUMBER_SPECIES_ARRAY = 11000` is the dimension of the species array.
```
  dim1 = [1]
  CALL h5LTmake_dataset_f(file_id, 'n_layer', 1, dim1, H5T_STD_I32LE, jmax, error)
  CALL h5LTmake_dataset_f(file_id, 'n_species', 1, dim1, H5T_STD_I32LE, NUMBER_SPECIES_ARRAY, error)
  dim1 = [jmax]
  CALL h5LTmake_dataset_f(file_id, 'temp', 1, dim1, H5T_IEEE_F64LE, temp, error)
  CALL h5LTmake_dataset_f(file_id, 'pres', 1, dim1, H5T_IEEE_F64LE, pres, error)
  CALL h5LTmake_dataset_f(file_id, 'rho' , 1, dim1, H5T_IEEE_F64LE, rho , error)
  dim2 = [NUMBER_SPECIES_ARRAY,jmax]
  CALL h5LTmake_dataset_f(file_id, 'ndens', 2, dim2, H5T_IEEE_F64LE, n, error)
```

### Indices for species in the number density array, `ndens`
The allocation of indices for species within the number density array, ndens, is currently broader than necessary. This will be optimized in the future.
| index | species   |
| ----- | --------- |
| 1     | total     |
| 10    | e-        |
| 100   | H         |
| 101   | H+        |
| 200   | He        |
| 201   | He+       |
| 202   | He++      |
| 300   | Li        |
| 301   | Li+       |
| 302   | Li++      |
| 303   | Li+++     |
| ...   | ...       |
| 9999  | Es99+     |
| 10001 | H2O       |
| 10002 | CO2       |
| 10003 | O3        |
| 10004 | N2O       |
| 10005 | CO        |
| 10006 | CH4       |
| 10007 | O2        |
| 10008 | NO        |
| 10009 | SO2       |
| 10010 | NO2       |
| 10011 | NH3       |
| 10012 | HNO3      |
| 10013 | OH        |
| 10014 | HF        |
| 10015 | HCl       |
| 10016 | HBr       |
| 10017 | HI        |
| 10018 | ClO       |
| 10019 | OCS       |
| 10020 | H2CO      |
| 10021 | HOCl      |
| 10022 | N2        |
| 10023 | HCN       |
| 10024 | CH3Cl     |
| 10025 | H2O2      |
| 10026 | C2H2      |
| 10027 | C2H6      |
| 10028 | PH3       |
| 10029 | COF2      |
| 10030 | SF6       |
| 10031 | H2S       |
| 10032 | HCOOH     |
| 10033 | HO2       |
| 10034 | OH+       |
| 10035 | ClONO2    |
| 10036 | NO+       |
| 10037 | HOBr      |
| 10038 | C2H4      |
| 10039 | CH3OH     |
| 10040 | CH3Br     |
| 10041 | CH3CN     |
| 10042 | CF4       |
| 10043 | C4H2      |
| 10044 | HC3N      |
| 10045 | H2        |
| 10046 | CS        |
| 10047 | SO3       |
| 10048 | C2N2      |
| 10049 | COCl2     |
| 10050 | VO        |
| 10051 | AlO       |
| 10052 | MgO       |
| 10053 | TiO       |
| 10054 | SiO       |
| 10055 | CaO       |
| 10056 | YO        |
| 10057 | NH        |
| 10058 | CH        |
| 10059 | SiH       |
| 10060 | SH        |
| 10061 | PH        |
| 10062 | SiH2      |
| 10063 | SiO2      |
| 10064 | MgH       |
| 10065 | NaH       |
| 10066 | AlH       |
| 10067 | CrH       |
| 10068 | CaH       |
| 10069 | BeH       |
| 10070 | TiH       |
| 10071 | FeH       |
| 10072 | LiH       |
| 10073 | ScH       |
| 10074 | SiH4      |
| 10075 | CH3F      |
| 10076 | AsH3      |
| 10077 | P2H2      |
| 10078 | PF3       |
| 10079 | CH3       |
| 10080 | PO        |
| 10081 | PN        |
| 10082 | KCl       |
| 10083 | NaCl      |
| 10084 | LiCl      |
| 10085 | CN        |
| 10086 | C2        |
| 10087 | CP        |
| 10088 | PS        |
| 10089 | NS        |
| 10090 | SiS       |
| 10091 | NaF       |
| 10092 | AlCl      |
| 10093 | AlF       |
| 10094 | KF        |
| 10095 | LiF       |
| 10096 | CaF       |
| 10097 | MgF       |
| 10098 | LiH+      |
| 10099 | H2+       |
| 10100 | HeH+      |
| 10101 | H3+       |
