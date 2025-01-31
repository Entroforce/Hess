# Hess

Hess is a platform for computational chemistry. This public version of
the platform provides functionality for empirical molecular docking.
It uses BFGS for local optimization and two different forms of Monte
Carlo method for global optimization. The Vinardo scoring function is
used for energy estimation.

## Prerequisities

* GNU/Linux operating system.
* `gcc`
* `g++`
* `make`
* `libz-dev`
* `libeigen3-dev`

## Build

* Activate the Intel compiler environment for `ifort`:
```
source /opt/intel/oneapi/setvars.sh
```
(skip this if you do not need Hess SE functionality). 

* Run `./setup.sh`, or `./setup-no-se.sh` if you do not need Hess SE functionality.

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
`--optimize arg`                global search algorithm variant (`mc` for Monte Carlo optimization or `mc_metropolis` for Monte Carlo with Metropolis criterion) <br/>
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

`./dist/Release/GNU-Linux/hess-empirical-docking -r protein_example.pdb -l ligand_example.pdb --autobox_ligand crystal_example.pdb --depth 500 --number_of_iterations 40 -o results.sdf`

## Preparation

### Usage

hess-preparation `<path to input pdb file>` `<path to output molfile>`

### Example

Go to the `hess-preparation` folder.

`./dist/Release/GNU-Linux/hess-preparation molecule_example.pdb prepared_molecule.mol`
