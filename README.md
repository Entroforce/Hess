# Hess

Hess is a toolkit for computational chemistry. This public version of
the toolkit provides functionality for empirical and semi-empirical
molecular docking.

## Prerequisities

* GNU/Linux operating system.
* `g++`
* `make`
* `libz-dev`
* `libeigen3-dev`
* Intel OneAPI toolkit with `icx` and `ifort` compilers (if you need Hess SE)

## Build

* Activate the Intel compiler environment for `ifort`:
```
source /opt/intel/oneapi/setvars.sh
```
(skip this if you do not need Hess SE). 

* Run `./setup.sh`, or `./setup-no-se.sh` if you do not need Hess SE.

(Indigo library will be downloaded in the process)

## Empirical docking (includes preparation)

### Usage

**Input:** (`pdb` or `mol` formats)<br/>
`-r [--receptor ] arg`          path to the receptor molecule<br/>
`-l [ --ligand ] arg`           path to the ligand molecule <br/>
`--autobox_ligand arg`          ligand to use for autobox <br/>

**Search space:** <br/>
`--center_x arg`                X coordinate of the center <br/>
`--center_x arg`                Y coordinate of the center <br/>
`--center_x arg`                Z coordinate of the center <br/>
`--size_x arg`                  size in the X dimension (Angstroms) <br/>
`--size_y arg`                  size in the Y dimension (Angstroms) <br/>
`--size_z arg`                  size in the Z dimension (Angstroms) <br/>

**Minimization options:** <br/>
`--optimize arg`                global search algorithm variant (`mc_metropolis` for Monte Carlo with Metropolis criterio or `mc` for Monte Carlo optimisation) <br/>
`--top arg`                     number of top scored poses stred in the output<br/>
`--depth arg`                   search depth (number of LBFGS local search runs) <br/>
`--granularity arg`             grid granularity (the smaller, the better the approximation) <br/>
`--score_only`                  score provided ligand pose <br/>
`--grid_deriv`                  use a grid with precalculated gradients <br/>

**Output:**  <br/>
`-o arg`                        output file name (`sdf`)<br/>

**Misc:** <br/>
`--seed arg`                    explicit random seed <br/>

### Example run

Go to the `hess-empirical-docking` folder.

`./dist/Release/GNU-Linux/hess-empirical-docking -r ../data/protein_example.pdb -l ../data/ligand_example.pdb --autobox_ligand ../data/crystal_example.pdb --number_of_iterations 16 -o results.sdf`

## Preparation

### Usage

hess-preparation `<path to input pdb file>` `<path to output molfile>`

### Example

Go to the `hess-preparation` folder.

`./dist/Release/GNU-Linux/hess-preparation ../data/molecule_example.pdb prepared_molecule.mol`

## SE program

SE program calculates the energy of a molecule. The program does not include the preparation of the molecule.

### Usage

For non-optimization run:

./hess-se-docking -p `<path to input pdb file>` -popt n -pout `<path to output molfile>`

For optimization run:

./hess-se-docking -p `<path to input pdb file>` -popt y -pout `<path to output molfile>`

### Example

Go to the `hess-se-docking` folder and run for 6MDC complex:

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_c_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_c_h_out.mol`

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_l_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_l_h_out.mol`

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_p_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_p_h_out.mol`

As a result, at each launch the energy of the molecule will be given out; to obtain the binding energy, the energy of the protein and ligand should be subtracted from the energy of the complex.

The same can be done for the 6S9B complex, but now with optimization enabled. Optimizing the protein and complex will take considerable time.

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_c_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_c_h_out.mol`

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_l_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_l_h_out.mol`

`./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_p_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_p_h_out.mol`
