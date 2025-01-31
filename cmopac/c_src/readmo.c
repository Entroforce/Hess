
#include <stdlib.h>
#include <string.h>

extern int nvar;
extern char keywrd[1000];
extern double tleft;

void setFpcByIdx(int *i, double *val);
void setMethodsByIdx(int *i, int *val);

extern double fpcref[20];

void getCoord(double **);
void getLabels(int **);
void getGeo(double **);
void getLoc(int **);
void getLopt(int **);
void getNa(int **);
void getNb(int **);
void getNc(int **);
void getXparam(double **);
void getAtmass(double **);
void getAms(double **);


extern void wrtkey_();

void readmo(int natoms_, double *coord_, int *labels_, const char* keywrd_, double *geo_int, int *loc_, int *lopt_, int *na_, int *nb_, int *nc_, int nvar_, double *xparam_){

    double *coord;
    getCoord(&coord);
    
    for(int i = 0; i < strlen(keywrd_); i++){
        keywrd[i] = keywrd_[i];
    }

    
    for(int i = 0; i < 10; i++){
        setFpcByIdx(&i, &(fpcref[2*i]));
    }
    
    double *geo;
    getGeo(&geo);

    int *lopt;
    getLopt(&lopt);
    for(int i = 0; i < natoms_ * 3; i++){
        geo[i] = geo_int[i];
        coord[i] = coord_[i];
        lopt[i] = lopt_[i];
    }
    if(nvar_ > 0){
        int *loc;
        getLoc(&loc);
        double *xparam;
        getXparam(&xparam);
        for(int i = 0; i < natoms_ * 3; i++){
            loc[2 * i] = loc_[2 * i];
            loc[2 * i + 1] = loc_[2 * i + 1];
            xparam[i] = xparam_[i];
        }
    }
    
    int *labels, *na, *nb, *nc;
    getLabels(&labels);
    getNa(&na);
    getNb(&nb);
    getNc(&nc);
    double *atmass;
    getAtmass(&atmass);
    double *ams;
    getAms(&ams);
    for(int i = 0; i < natoms_; i++){
        labels[i] = labels_[i];
        atmass[i] = ams[labels[i]-1];
        na[i] = na_[i];
        nb[i] = nb_[i];
        nc[i] = nc_[i];
    }
    nvar = nvar_;

    tleft = 172800.0;
    
    int pm7_one = 1;
    int pm7_idx = 13;
    if(strstr(keywrd_, "PM7") == NULL){
        pm7_idx = 2;
    }
    setMethodsByIdx(&pm7_idx, &pm7_one);
    
    wrtkey_();
}