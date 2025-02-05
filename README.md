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

### With Hess SE

```
source /opt/intel/oneapi/setvars.sh

./setup.sh
```

Indigo library will be downloaded in the process.

### Without Hess SE

```
./setup-no-se.sh
```

Indigo library will be downloaded in the process.

## Empirical docking (includes preparation)

### Usage

```
hess-empirical-docking [options]
```

#### Input options 

PDB and Molfile formats are accepted.

- `-r, --receptor arg`:  path to the receptor molecule
- `-l, --ligand arg`: path to the ligand molecule
- `--autobox_ligand arg`: ligand to use for autobox calculation

#### Box options

Units are Angstroms.

- `--center_x arg`:                X coordinate of the center
- `--center_x arg`:                Y coordinate of the center
- `--center_x arg`:                Z coordinate of the center 
- `--size_x arg`:                  size in the X dimension
- `--size_y arg`:                  size in the Y dimension
- `--size_z arg`:                  size in the Z dimension

#### Minimization algorithm options

- `--optimize arg`:                global search algorithm variant (`mc_metropolis` for Monte Carlo with Metropolis criterion or `mc` for Monte Carlo optimisation) 
- `--top arg`:                     number of top scored poses stred in the output
- `--depth arg`:                   search depth (number of LBFGS local search runs)
- `--granularity arg`:             grid granularity (the smaller, the better the approximation)
- `--score_only`:                  score provided ligand pose
- `--grid_deriv`:                  use a grid with precalculated gradients

#### Output

- `-o arg`:                        output file name (`sdf`)

#### Other

- `--seed arg`:                    explicit random seed

### Example run

Go to the `hess-empirical-docking` folder.

```
./dist/Release/GNU-Linux/hess-empirical-docking -r ../data/protein_example.pdb -l ../data/ligand_example.pdb --autobox_ligand ../data/crystal_example.pdb --number_of_iterations 16 -o results.sdf
```

## Preparation

### Usage

```
hess-preparation <path to input pdb file> <path to output molfile>
```

### Example

Go to the `hess-preparation` folder.

`./dist/Release/GNU-Linux/hess-preparation ../data/molecule_example.pdb prepared_molecule.mol`

## SE docking

### Usage

For non-optimization run:

```
./hess-se-docking -p <path to input pdb file> -popt n -pout <path to output molfile>
```

For optimization run:

```
./hess-se-docking -p <path to input pdb file> -popt y -pout <path to output molfile>
```

### Example

#### Without optimization

Go to the `hess-se-docking` folder and run the energy calculation for 6MDC complex:

```
./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_c_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_c_h_out.mol
```

And then for the ligand and the protein separately:

```
./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_l_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_l_h_out.mol

./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6mdc/6mdc_p_h.pdb -popt n -pout ../data/test-se/6mdc/6mdc_p_h_out.mol
```

After each run, the energy will be printed out; to obtain the binding energy, the energy of the protein and ligand should be subtracted from the energy of the complex.

#### With optimization

Here is a similar example for the 6S9B complex, but now with optimization enabled. Optimizing the protein and complex will take considerable time.

```
./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_c_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_c_h_out.mol

./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_l_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_l_h_out.mol

./dist/Release/GNU-Linux/hess-se-docking -p ../data/test-se/6s9b/6s9b_p_h.pdb -popt y -pout ../data/test-se/6s9b/6s9b_p_h_out.mol
```
