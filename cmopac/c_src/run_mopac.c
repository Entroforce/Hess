#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "mkl.h"

extern void dhc_(double *, double *, double *, double *, double *, int*, int*, int*, int*, double*, int*, int*);
extern void fock2_(double *f, double *ptot, double *p, double *w, double *wj, double *wk, int *numat, int *nfirst, int *nlast, int *mode, int *isgood);
extern void deriv_(double *geo, double *grad_, int *isgood);
extern void dcart_(double *geo, double *grad_, int* isgood);
extern void picopt_(int *minus_one);
extern void post_scf_corrections_(double * corrects, int *scf_false);
extern void setup_mopac_arrays_(int *n, int *mode);
extern void delete_mozyme_arrays_();
extern void summary_(char *txt, int *ntxt);
extern void switch_();
extern void moldat_(int *zero);
extern void calpar_();
extern double c_triple_bond_c_();
extern void lbfgs_(double *xparam_, double *funct);
extern void fbx_();
extern void fordd_();
extern void geochk_();
extern void set_up_mozyme_arrays_();
extern void compfg(double *xparam_, int int_bool, double *escf, int fulscf, double *grad_, int lgrad_);
extern void chrge_(double *p, double *q);
extern void iter_for_mozyme_(double * elect_, int *isgood);
extern void gmetry_(double *geo, double *coord);
extern void readmo(int natoms, double *coord, int* labels, const char* keywrd_, double *geo_int, int *loc, int *lopt, int *na, int *nb, int *nc, int nvar_, double *xparam_);

//void getCommonArraysC(int *);

// common_arrays_c
void getNfirst(int **);
void getNlast(int **);
void getNat(int **);
void getXparam(double **);
void getGrad(double **);
void getNw(int **);
void getP(double **);
void getTxtatm(char **);
void getGeo(double **);
void getQ(double **);


// molkst_c
extern double gnorm, escf, time0, atheat, stress, e_disp,
        e_hb, e_hh, param_constant, trunc_1, trunc_2, elect, enuclr,
        sum_dihed, press[3];

void setMethodsByIdx(int *i, int *val);
void setMM_correctionsByIdx(int *i, int *val);

extern int natoms, numat, nvar, numcal, nscf, 
           iflepo, iscf, last, moperr, maxatoms,
        isok, na1, mozyme, step_num, no_pka, lxfac, in_house_only,
        sparkle, num_threads, uhf, nbeta, time_start[8];

extern char keywrd[1000], errtxt[200], allkey[1000];

// parameters_c
extern int ios[107], iop[107], iod[107];
extern double eisol[107], eheat[107], zs[107], tore[107];

extern char elemnt[107][2];

extern int iseps;
extern int useps;
extern int lpka;
extern double solv_energy;

extern double rxn_coord;
extern char state_irred_rep[4];

extern char density_fn[241];

extern int nres;

extern int nmos;

extern int lgpu;


void end_function() {
    int zero = 0;
    int one = 1;
    char *space = " ";
    
    for(int i = 0; i < 17; i++){
        setMethodsByIdx(&i, &zero);
    }
    for(int i = 0; i < 1000; i++){
        allkey[i] = ' ';
        keywrd[i] = ' ';
    }

    nbeta = 0;
    uhf = 0;
    
    setup_mopac_arrays_(&zero, &zero);
    delete_mozyme_arrays_();
    summary_(space, &one);
    // printf("==MOPAC DONE==\n");
    return;
}

double run_mopac(double *coord, int *labels, char **txtatm_, int natoms_, const char *keywrd_, 
        double *geo_int, int *loc, int *lopt, int *na, int *nb, 
        int *nc, int nvar_, double *xparam_, double *coords, double *q_, double *escf_, double *grad_) {
    moperr = 0;
    printf("run with keywords: %s\n", keywrd_);
    for (int i = 0; i < 107; i++) {
        tore[i] = ios[i] + iop[i] + iod[i];
    }
    
    fbx_();    // Factorials and Pascal's triangle (pure constants)
    fordd_();  // More constants, for use by MNDO-d
    param_constant = 1.0;
    
    char filename[241] = "mol.den";
    
    memcpy(density_fn, filename, 241);
    
    lgpu = 0;

    if (param_constant < -1e05) {
        return 0.0;
    }

    trunc_1 = 7.0;   // Beyond 7.0 Angstroms, use exact point-charge
    trunc_2 = 0.22;  // Multiplier in Gaussian: exp(-trunc_2*(trunc_1 - Rab)^2)

    natoms = natoms_;

    if (natoms == 0 || moperr) {
        return 0.0;
    }

    isok = 1;

    char *errtxt_ = "Job stopped by operator";
    strcpy(errtxt, errtxt_);

    for (int i = 0; i < 8; i++) {
        time_start[i] = 0;
    }

    if (moperr) {
        end_function();
        return 0.0;
    }


    int one = 1;
    setup_mopac_arrays_(&natoms, &one);
    maxatoms = natoms;

    in_house_only = 0;
    numcal = 1;
    step_num = 1;
    moperr = 0;
    escf = 0.0;
    gnorm = 0.0;
    press[0] = 0.0;
    press[1] = 0.0;
    press[2] = 0.0;
    e_disp = 0.0;
    e_hb = 0.0;
    e_hh = 0.0;
    solv_energy = 0.0;
    nres = 0.0;
    nscf = 0.0;
    nmos = 0.0;
    na1 = 0.0;
    lpka = 0;
    stress = 0.0;
    // no_pka = 0;
    time0 = 0.0;
    
    //zero 
    int zero = 0;

    for (int i = 0; i < 20; i++) {
        setMM_correctionsByIdx(&i, &zero);
    }
    char *space = " ";
    strcpy(state_irred_rep, space);
    
    char *txtatm;
    getTxtatm(&txtatm);
    
    for(int i = 0; i < natoms_; i++){
        for(int j = 0; j < 21; j++){
            txtatm[26 * i + j] = txtatm_[i][j];
        }
    }

    readmo(natoms, coord, labels, keywrd_, geo_int, loc, lopt, na, nb, nc, nvar_, xparam_);

    if (moperr) {
        end_function();
    }

    if (numcal == 1) {
        num_threads = 8;
        printf("threads: %d\n", num_threads);
        mkl_set_num_threads(num_threads);
    }

    lxfac = 0;
    int thirteen = 13;
    setMethodsByIdx(&thirteen, &one);
    switch_();

    sparkle = 0;

    for (int i = 56; i < 71; i++) {
        if (zs[i] < 0.1) {
            tore[i] = 3.0;
        }
    }

    moldat_(&one);
    calpar_();

    useps = 0;
    iseps = strstr(keywrd, "EPS=") || strstr(keywrd, "PKA");
    mozyme = strstr(keywrd, "MOZYME") != NULL;

    int two = 2;
    setup_mopac_arrays_(&one, &two);
    //cosmo_c_mp_iseps_ = strstr(keywrd, "EPS=") != NULL;
    
    int *nw;
    getNw(&nw);
    if (nw) {
        // BAD!!! TODO!
        //free(nw);
    }
    nw = (int *)malloc(numat * sizeof(int));

    
    int *nfirst, *nlast, *nat;
    getNat(&nat);
    getNfirst(&nfirst);
    getNlast(&nlast);
    
    double *xparam, *grad, *p;
    getXparam(&xparam);
    getGrad(&grad);
    getP(&p);
    
    
    int l = 1;
    for (int i = 0; i < numat; i++) {
        nw[i] = l;
        l = l + ((nlast[i] - nfirst[i] + 1) * (nlast[i] - nfirst[i] + 2)) / 2;
    }

    // calculate atomic energy

    double eheat_sum = 0;
    double eat = 0;
    for (int i = 0; i < numat; i++) {
        eheat_sum += eheat[nat[i] - 1];
        eat += eisol[nat[i] - 1];
    }
    atheat = eheat_sum;
    atheat -= eat * 23.060529;
    atheat += c_triple_bond_c_();
    rxn_coord = 1.e9;
    int computed_charge;
    if (mozyme) {
        geochk_(&computed_charge);
        if(strstr(keywrd, "CHARGES")){
            end_function();
            return computed_charge;
        }
        set_up_mozyme_arrays_();
        int minus_one = -1;
        picopt_(&minus_one);
    }

    double some_value_for_mozyme_reload = 0.0;
    int isgood = 1;
    int some_zero = 0;
    iter_for_mozyme_(&some_value_for_mozyme_reload, &isgood);
    
    dhc_(NULL, NULL, NULL, NULL, NULL, &some_zero, &some_zero, &some_zero, &some_zero, &some_value_for_mozyme_reload, &some_zero, &isgood);
    dcart_(NULL,NULL,&isgood);

    if (strstr(keywrd, "1SCF") || nvar == 0) {
        iflepo = 1;
        iscf = 1;
        last = 1;
        for (int i = 0; i < nvar; i++) {
            grad[i] = 0.0;
        }
        numcal = 1;
        int scf_true = 1;
        int scf_false = 0;
        int deriv_one = 1;
        deriv_(NULL, NULL, &deriv_one);
        compfg(xparam, 1, &escf, 1, grad, scf_false);
    } 
    else if(strstr(keywrd, "LBFGS")){
        int deriv_one = 1;
        deriv_(NULL, NULL, &deriv_one);
        lbfgs_(xparam, &escf);
    }
    else{
        // printf("GEOMETRY OPTIMIZATION USING EF\n");
        //ef_(xparam, &molkst_c_mp_escf_);
    }
    last = 1;
    double corrects;
    int post_scf_false = 0;
    post_scf_corrections_(&corrects, &post_scf_false);
    
    double *geo, *q;
    getGeo(&geo);
    getQ(&q);
    
    gmetry_(geo, coords);
//    if (q_ != NULL && q!= NULL) {
//        chrge_(p, q_);
//    }

    
    if(grad_ != NULL){
        numat = natoms_;
        int good = 0;
        int testgood = 1;
        fock2_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &testgood);
        dcart_(coords, grad_, &good);
    }

    printf("coords:\n");
    for(int i = 0; i < natoms_; i++){
        printf("%c%c %lf %lf %lf\n", elemnt[nat[i]-1][0],elemnt[nat[i]-1][1], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
    }
    
    double ev_1_in_kcal_mol = 23.060529;
    *escf_ = escf/ev_1_in_kcal_mol;
    end_function();
    printf("Energy terms: \nElectronic energy: %lf eV \nNuclear-nuclear repulsion energy: %lf eV \nCorrects: %lf eV \nDihed: %lf eV \nSolv_energy: %lf eV \n", elect, enuclr, corrects/ev_1_in_kcal_mol, sum_dihed/ev_1_in_kcal_mol, solv_energy);          
    return elect + enuclr + corrects/ev_1_in_kcal_mol + sum_dihed/ev_1_in_kcal_mol + solv_energy;
}