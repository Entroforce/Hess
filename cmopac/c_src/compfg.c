#include <string.h>
#include <math.h>

void getNat(int **);
void getLoc(int **);
void getGeo(double **);
void getNa(int **);
void getNb(int **);
void getNc(int **);
void getCoord(double **);
void getXparef(double **);
void getTvec(double **);
void getLabels(int **);
void getF(double **);

void getMethodsByIdx(int *, int*);

extern int numat, norbs, nelecs, nclose, nopen, numcal, ndep, nvar, moperr, mozyme, id,
        iseps, useps, noeps, nnhco, nhco[4000][4];

extern double fract, elect, enuclr, emin, atheat, pressure, density, stress, hpress, sum_dihed,
                nsp2_corr, si_o_h_corr, solv_energy, htype;

extern char keywrd[1000];

void getPartf(double **);

void linear_cosmo_mp_ini_linear_cosmo_();
void linear_cosmo_mp_coscavz_();

extern void symtry_();
extern void hcore_();
extern void gmetry_(double *geo, double *coord);

extern double dot_(double *x,double *y,int *n);
extern double volume_(double *vecs, int *ndim);

extern double nsp2_correction_();
extern double si_o_h_correction_();

extern void dihed_(double *xyz, int *i, int *j, int *k, int *l, double *angle);
extern void post_scf_corrections_(double *correction, int *l_grad);

//void iter(double *ee, int fulscf, int rand_bool);

extern void deriv_(double *geo, double *grad, int *isgood);


extern void hcore_for_mozyme_();
extern void iter_for_mozyme_(double * elect, int * isgood);
extern double helecz_();
extern void buildf_(double *f, double *partf, int *mode);
extern void cosini_(int *a);
extern void mopend_(char *str);
extern void linear_cosmo_mp_ini_linear_cosmo_();
extern void linear_cosmo_mp_coscavz_(double *coord, int* nat);
extern void coscav_();
extern void mkbmat_();

int compfg_debug;
int compfg_print;
int compfg_large;
int compfg_usedci;
int compfg_force;
int compfg_times;
int compfg_aider;
int compfg_icalcn = 0;
int compfg_dh;
int compfg_dl_locate_ts;

void compfg(double *xparam, int int_bool, double * escf, int fulscf, double *grad, int lgrad){
    
    double angle, atheat_store, sum;
    compfg_icalcn = 0;
    if(compfg_icalcn != numcal){
        moperr = 0;
        compfg_icalcn = numcal;
        hpress = 0.0;
        nsp2_corr = 0.0;
        si_o_h_corr = 0.0;
        sum_dihed = 0.0;
        
        if(iseps){
            iseps = 1;
            noeps = 1;
            //cosmo_c_mp_useps_ = 1;
            int cos_true = 1;
            cosini_(&cos_true);
            if(moperr){
                return;
            }
            mozyme = strstr(keywrd, "MOZ") || strstr(keywrd, "LOCATE-TS") || strstr(keywrd, "RAPID") || strstr(keywrd, "REFINE-TS");
            if(mozyme && numat == 1){
                char *mopend_str = "MOZYME cannot be used for systems composed of only one atom!";
                mopend_(mopend_str);
                return;
            }
            if(mozyme){
                linear_cosmo_mp_ini_linear_cosmo_();
            }
        }
        
        
        compfg_aider = 0;
        compfg_times = 0;
        compfg_usedci = (nclose != nopen) &&
                 (fabs(fract - 2.0) > 1e-20) && 
                 (fract > 1e-20) || 0;

        compfg_force = 0;
        compfg_large = 0;
        compfg_debug = 0;
        compfg_print = 0;
        compfg_dl_locate_ts = 0;
        
        int mp_methods;
        int thirteen = 13;
        getMethodsByIdx(&thirteen, &mp_methods);
        compfg_dh = strstr(keywrd, " PM6-D") || strstr(keywrd, " PM6-H") || (mp_methods == 1);
        emin = 0.0;
        double *xparef;
        getXparef(&xparef);
        for(int i = 0; i < nvar; i++){
            xparef[i] = xparam[i];
        }
    }


    //set up coordinates for current calculation

    int *loc;
    getLoc(&loc);
    double *geo;
    getGeo(&geo);
    double *coord;
    getCoord(&coord);
    
    for (int i = 0; i < nvar; i++){
        geo[3 * (loc[2*i] - 1) + (loc[2*i+1] - 1)] = xparam[i];
    }
    if(ndep != 0){
        symtry_();
    }

    gmetry_(geo, coord);

    if(moperr) {
        return;
    }
    
    if(iseps){
        if(mozyme){
            //error here
            int *nat;
            getNat(&nat);
            linear_cosmo_mp_coscavz_(coord, nat);
        }
        else {
            coscav_();
            mkbmat_();
        }
        if(moperr) {
            return;
        }
        if(noeps){
            useps = 0;
        }
    }
    
    if(mozyme){
        if(iseps)
            useps = 1;
        if(compfg_dl_locate_ts || int_bool)
            hcore_for_mozyme_();
    }
    else if(int_bool){
        hcore_();       
    }

    if(moperr){
        return;
    }

    atheat_store = atheat;

    if((norbs > 0) && (nelecs > 0)){
        hpress = 0.0;
        if(fabs(pressure) > 1e-4){
            int three = 3;
            double *tvec;
            getTvec(&tvec);
            if(id == 1) {
                hpress = -pressure * sqrt(dot_(tvec, tvec, &three));
            }
            else if(id == 3){
                hpress = -pressure * volume_(tvec, &three);
            }
            atheat += hpress;
        }

        if(useps && !mozyme){
            atheat += solv_energy * 23.060529;
        }

        if(0){
            nsp2_corr = nsp2_correction_();
            atheat += nsp2_corr;
        }

//         if(molkst_c_mp_methods_[13] == 1 && molkst_c_mp_si_o_h_present_) {
//             molkst_c_mp_si_o_h_corr_ = si_o_h_correction_();
//             molkst_c_mp_atheat_ += molkst_c_mp_si_o_h_corr_;
//         }

        sum_dihed = 0.0;

        for(int i = 0; i < nnhco; i++){
            dihed_(coord, &nhco[i][0], &nhco[i][1], &nhco[i][2], &nhco[i][3], &angle);
            sum_dihed += htype * pow(sin(angle),2);
        }

        atheat += sum_dihed;
        stress = 0.0;

        int thirteen = 13;
        int mp_methods;
        getMethodsByIdx(&thirteen, &mp_methods);
        
        if(compfg_dh && (mp_methods == 1)) {
            int post_scf_false = 0;
            post_scf_corrections_(&sum, &post_scf_false);
            if(moperr){
                return;
            }
            atheat += sum;
        }

        if(int_bool){
            if(mozyme){
                int isgood = 0;
                iter_for_mozyme_(&elect, &isgood);
            }
            else {
                //iter(&molkst_c_mp_elect_, fulscf, 1);
                //bool iter_true = true;
                //iter_(&molkst_c_mp_elect_, &fulscf, &iter_true);
            }
            if(moperr) {
                return;
            }
            if(noeps){
                noeps = 0;
                useps = 1;
                if(!mozyme){
                    hcore_();
                    //iter(&molkst_c_mp_elect_, fulscf, 1);
                    //iter_(&molkst_c_mp_elect_, &fulscf, &iter_true);
                }
            }
        } else {
            if(mozyme){
                int zero = 0;
                double *partf;
                getPartf(&partf);
                double *f;
                getF(&f);
                buildf_(f, partf, &zero);
                elect = helecz_();
            }
        }

        stress *= density;
        atheat += stress;
        if(moperr){
            return;
        }
    }
    else {
        elect = 0.0;
    }
    //23.060529 -- some physical constant, fpc9
    *escf = (elect + enuclr) * 23.060529 + atheat;
    if(useps && mozyme){
        *escf = *escf + solv_energy * 23.060529;
    }
    
    int thirteen = 13;
    int mp_methods;
    getMethodsByIdx(&thirteen, &mp_methods);
    
    if(compfg_dh && (mp_methods != 1)){
        int post_scf_false = 0;
        post_scf_corrections_(&sum, &post_scf_false);
        if(moperr){
            return;
        }
        *escf = *escf + sum;
    }
    atheat = atheat_store;
    if(*escf < emin || emin == 0.0){
        emin = *escf;
    }

    if(lgrad){
        if(nelecs > 0){
            int deriv_zero = 0;
            deriv_(geo, grad, &deriv_zero);
        }
        if(moperr){
            return;
        }
    }

    return;
}
