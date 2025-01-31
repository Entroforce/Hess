#include <stdio.h>
#include <stdbool.h>
typedef struct cnfg {
    char *args;
    FILE *input;
    char *keywrd;
    int natoms;
    int maxatoms;
    int numcal;
    int step_num;
    int **lopt; //optimization flags
    bool moperr;
    char *txtatm1;
    double **break_coords;
    // double escf;
    // double gnorm;
    // double press;
    // double E_disp;
    // double E_hb;
    // double E_hh;
    // double solv_energy;
    // int nres;
    // int nscf;
    // int nmos;
    // int na1;
    // bool lpka;
    // double stress;
    // int no_pKa;
    // time0 = second(1)
    // bool MM_corrections;
    // char *state_Irred_Rep;

} Config;
