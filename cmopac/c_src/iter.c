#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

extern int molkst_c_mp_n2elec_;

extern double *common_arrays_c_mp_eigs_;
extern double *common_arrays_c_mp_p_;
extern double *common_arrays_c_mp_pa_;
extern double *common_arrays_c_mp_pb_;
extern double *common_arrays_c_mp_cb_; //[norbs][norbs]
extern double *common_arrays_c_mp_h_;
extern double *common_arrays_c_mp_c_; //[norbs][norbs]
extern int *common_arrays_c_mp_nat_;
extern int *common_arrays_c_mp_nfirst_;
extern int *common_arrays_c_mp_nlast_;
extern double *common_arrays_c_mp_eigb_;
extern double *common_arrays_c_mp_pdiag_;
extern double *common_arrays_c_mp_f_;
extern double *common_arrays_c_mp_w_;
extern double *common_arrays_c_mp_wk_;
extern double *common_arrays_c_mp_fb_;

extern double *funcon_c_mp_fpc_;

extern double *iter_c_mp_pold_;
extern double *iter_c_mp_pold2_;
extern double *iter_c_mp_pbold_;
extern double *iter_c_mp_pbold2_;
extern double *iter_c_mp_pold3_;
extern double *iter_c_mp_pbold3_;
extern double *iter_c_mp_vec_ai_;  //[norbs][norbs]
extern double *iter_c_mp_vec_bi_;  //[norbs][norbs]
extern double *iter_c_mp_fock_ai_; //[norbs][norbs]
extern double *iter_c_mp_fock_bi_; //[norbs][norbs]
extern double *iter_c_mp_p_ai_;    //[norbs][norbs]
extern double *iter_c_mp_p_bi_;    //[norbs][norbs]
extern double *iter_c_mp_h_ai_;
extern double *iter_c_mp_h_bi_;
extern double *iter_c_mp_vecl_ai_;
extern double *iter_c_mp_vecl_bi_;

extern double funcon_c_mp_fpc_9_;

extern int maps_c_mp_latom_;

extern int chanel_c_mp_iw_;
extern int chanel_c_mp_ifiles_1_;

extern int molkst_c_mp_numat_;
extern int molkst_c_mp_norbs_;
extern int molkst_c_mp_nalpha_;
extern int molkst_c_mp_nbeta_;
extern int molkst_c_mp_uhf_;
extern int molkst_c_mp_nclose_;
extern int molkst_c_mp_nopen_;
extern double molkst_c_mp_fract_;
extern int molkst_c_mp_numcal_;
extern int molkst_c_mp_mpack_;
extern int molkst_c_mp_iflepo_;
extern int molkst_c_mp_iscf_;
extern double molkst_c_mp_enuclr_;
extern char molkst_c_mp_keywrd_[1000];
extern double molkst_c_mp_gnorm_;
extern int molkst_c_mp_moperr_;
extern int molkst_c_mp_last_;
extern int molkst_c_mp_nscf_;
extern double molkst_c_mp_emin_;
extern int molkst_c_mp_limscf_;
extern double molkst_c_mp_atheat_;
extern int molkst_c_mp_is_param_;
extern int molkst_c_mp_id_;
extern char molkst_c_mp_line_[1000];
extern int molkst_c_mp_lxfac_;
extern int molkst_c_mp_nalpha_open_;
extern int molkst_c_mp_nbeta_open_;

extern void delete_iter_arrays_();
extern void dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
extern void vecprt_(double *A, int *numm);

extern void fock2_(double *f, double *ptot, double *p, double *w, double *wj, double *wk, int *numat, int *nfirst, int *nlast, int *mode);
void fock2(double *f, double *ptot, double *p, double *w, double *wj, double *wk,
           int numat, int *nfirst, int *nlast, int mode, int *ifact, int *i1fact,
           double **ptot2, int *iter_icalcn, int *ione, int *jindex, int *lid);

extern void mopend_(char *txt);

extern double helect_(int *N, double *P, double *H, double *F);
double helect(int N, double *P, double *H, double *F);


extern double capcor_(int *NAT, int *NFIRST, int *NLAST, double *P, double *H);

extern void interp_(int *np, int *nq, int *mode, double *e, double *fp, double *cp, double *theta, double *vec_interp, double *fock_interp, double *p_interp, double *h_interp, double *vecl, double *eold_ref);
extern void pulay_(double *F, double *P, int *N, double *FPPF, double *FOCK, double *EMAT, int *LFOCK, int *NFOCK, int *MSIZE, int *START, double *iter_pl);

extern void eigenvectors_lapack_(double *eigenvecs, double *xmat, double *eigvals, int *ndim);

extern void diag_for_gpu_(double *fao, double *vector, int *nocc, double *eig, int *norbs, int *mpack);

extern void swap_(double *C, int *N, int *MDIM, int *NOCC, int *iter_ifill);

extern void density_for_gpu_(double *c, double *fract, int *ndubl, int *nsingl, double *occ, int *mpack, int *norbs, int *mode, double *pp, int *iopc);

extern void cnvg_(double *PNEW, double *P, double *P1, int *NITER, double *iter_pl);

extern void densit_(double *c, int *mdim, int *norbs, int *nocc, double *occ, int *nfract, double *fract, double *p, int *mode);

extern void phase_lock_(double *vecs, int *n);

extern double meci_();

int iter_icalcn = 0;
int iter_debug = 0;
int iter_prtfok = 0;
double iter_titer0 = 0.0;
int iter_prteig = 0;
int iter_prtden = 0;
int iter_prt1el = 0;
double iter_ten = 10.0;
double iter_tenold = 10.0;
double iter_plb = 0.0;
double iter_scorr = 0.0;
char *iter_abprt[3][6] = {"     ", "ALPHA", " BETA"};

int iter_minprt;
int iter_newdg;
double iter_scfcrt;
int iter_prtpl;
int iter_prtvec;
double iter_pl;
double iter_bshift;
double iter_pltest;
int iter_itrmax;
int iter_na2el;
int iter_na1el;
int iter_nb1el;
int iter_ifill;
int iter_camkin;
int iter_ci;
int iter_okpuly;
int iter_oknewd;
int iter_times;
int iter_force;
int iter_allcon;
double iter_trans;
int iter_halfe;
double iter_w1, iter_w2;
double iter_random;
int iter_gs;
double iter_shift;
double iter_shiftb;
double iter_shfmax;
int iter_capps;
int iter_incitr;
int iter_irrr;
double iter_plchek;
int iter_timitr;
double iter_shfto;
double iter_shftbo;
int iter_jalp;
int iter_ialp;
int iter_jbet;
int iter_ibet;
double iter_enrgy;
int iter_iopc_calcp;

void iter(double *ee, int fulscf, int rand_bool)
{

    int *ifact;
    int *i1fact;
    double **ptot2;
    ptot2 = (double **)malloc(81 * sizeof(double *));
    for (int i = 0; i < 81; i++)
    {
        ptot2[i] = (double *)malloc(molkst_c_mp_numat_ * sizeof(double));
    }
    ifact = (int *)malloc((3 + molkst_c_mp_norbs_) * sizeof(int));
    i1fact = (int *)malloc((3 + molkst_c_mp_norbs_) * sizeof(int));
    int dd = 0;
    double selcon;
    int l, ihomo, ihomob, iemin,
        iemax, iredy, niter, modea, modeb;

    double q[molkst_c_mp_numat_];

    double escf0[10];

    double eold, diff, titer, escf,
        sellim, sum, eold_alpha, eold_beta;
    double theta[molkst_c_mp_norbs_];

    int frst, bfrst, ready, glow,
        makea, makeb, getout, l_param;

    int icalcn_fock2;
    int ione_fock2;
    int jindex_fock2[256];
    int lid_fock2;

    icalcn_fock2 = 0;

    iter_ifill = 0;
    ihomo = molkst_c_mp_nclose_ + molkst_c_mp_nalpha_ > 1 ? molkst_c_mp_nclose_ + molkst_c_mp_nalpha_ : 1;
    ihomob = molkst_c_mp_nclose_ + molkst_c_mp_nbeta_ > 1 ? molkst_c_mp_nclose_ + molkst_c_mp_nbeta_ : 1;
    eold = 100.0;
    ready = 0;
    diff = 0.0;
    for (int i = 0; i < 10; i++)
    {
        escf0[i] = 0.0;
    }

    if (iter_icalcn != molkst_c_mp_numcal_)
    {
        delete_iter_arrays_();
        l_param = 1;
        iter_enrgy = 23.060529;
        glow = 0;
        iter_irrr = 5;
        iter_shift = 0.0;
        iter_icalcn = molkst_c_mp_numcal_;
        iter_shfmax = 20.0;

        // iter_debug key-words worked out
        iter_debug = 0;
        iter_minprt = (maps_c_mp_latom_ == 0);
        iter_prteig = 0;
        iter_prtpl = strstr(molkst_c_mp_keywrd_, "PL") || strstr(molkst_c_mp_keywrd_, "PLS");
        iter_prt1el = 0;
        iter_prtden = 0;
        iter_prtfok = 0;
        iter_prtvec = 0;
        iter_debug = 0;

        // initialize some logicals and constants

        iter_newdg = 0;
        iter_camkin = 0;
        iter_plchek = 0.005;
        iter_pl = 1.0;
        iter_plb = 0.0;
        iter_bshift = -80.0;
        iter_shift = 1.0;
        iter_shfto = 0.0;
        iter_shftbo = 0.0;
        iter_itrmax = 2000;
        iter_na2el = molkst_c_mp_nclose_;
        iter_na1el = molkst_c_mp_nalpha_ + molkst_c_mp_nopen_;
        iter_nb1el = molkst_c_mp_nbeta_ + molkst_c_mp_nopen_;

        // use key-words to assign various constants

        if (strstr(molkst_c_mp_keywrd_, "SHIFT"))
        {
            iter_bshift = -10.0;
        }
        if (fabs(iter_bshift) > 1e-20)
        {
            iter_ten = iter_bshift;
        }
        if (strstr(molkst_c_mp_keywrd_, "ITRY"))
        {
            iter_itrmax = 2000;
        }
        iter_ci = 0;
        iter_okpuly = 0;
        iter_oknewd = fabs(iter_bshift) < 0.001;
        if (iter_camkin && (fabs(iter_bshift) > 1e-5))
        {
            iter_bshift = 4.44;
        }
        iter_times = 0;
        iter_timitr = iter_times;
        iter_force = 0;
        iter_gs = 1;
        iter_allcon = iter_okpuly || iter_camkin;

        int nat_counter = 0;
        for (int i = 0; i < molkst_c_mp_numat_; i++)
        {
            if (common_arrays_c_mp_nat_[i] == 102)
            {
                nat_counter += 1;
            }
        }

        iter_capps = (nat_counter > 0);
        molkst_c_mp_iscf_ = 1;
        iter_trans = 0.200;
        if (!molkst_c_mp_is_param_)
        {
            for (int i = 0; i < molkst_c_mp_mpack_; i++)
            {
                common_arrays_c_mp_p_[i] = 0.0;
                common_arrays_c_mp_pa_[i] = 0.0;
                common_arrays_c_mp_pb_[i] = 0.0;
            }
            iter_w1 = iter_na1el / (iter_na1el + 1e-6 + iter_nb1el);
            iter_w2 = 1.0 - iter_w1;
            if (iter_w1 < 1e-6)
            {
                iter_w1 = 0.5;
            }
            if (iter_w2 < 1e-6)
            {
                iter_w2 = 0.5;
            }
            iter_random = 1.0;
            glow = glow || (molkst_c_mp_gnorm_ < 2.0) && (molkst_c_mp_gnorm_ > 1e-9);
            if (!glow && molkst_c_mp_uhf_ && (iter_na1el == iter_nb1el))
            {
                iter_random = 1.1;
            }
            for (int i = 0; i < molkst_c_mp_norbs_; i++)
            {
                int j = ((i + 1) * (i + 2)) / 2 - 1;
                common_arrays_c_mp_p_[j] = common_arrays_c_mp_pdiag_[i];
                common_arrays_c_mp_pa_[j] = common_arrays_c_mp_p_[j] * iter_w1 * iter_random;
                iter_random = 1.0 / iter_random;
                common_arrays_c_mp_pb_[j] = common_arrays_c_mp_p_[j] * iter_w2 * iter_random;
            }
            if (molkst_c_mp_uhf_)
            {
                for (int i = 0; i < molkst_c_mp_norbs_; i++)
                {
                    iter_random = 1.0 / iter_random;
                    int j = ((i + 1) * (i + 2)) / 2 - 1;
                    common_arrays_c_mp_pb_[j] = common_arrays_c_mp_p_[j] * iter_w2 * iter_random;
                }
            }
        }

        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            iter_c_mp_pold_[i] = common_arrays_c_mp_pa_[i];
        }
        if (molkst_c_mp_uhf_)
        {
            for (int i = 0; i < molkst_c_mp_mpack_; i++)
            {
                iter_c_mp_pbold_[i] = common_arrays_c_mp_pb_[i];
            }
        }
        for (int i = 0; i < molkst_c_mp_norbs_; i++)
        {
            int j = ((i + 1) * (i + 2)) / 2 - 1;
            iter_c_mp_pold2_[i] = iter_c_mp_pold_[j];
        }
        iter_halfe = ((molkst_c_mp_nopen_ != molkst_c_mp_nclose_) && (fabs(molkst_c_mp_fract_ - 2.0) > 1e-20) && (fabs(molkst_c_mp_fract_) > 1e-20));
        if (iter_halfe)
        {
            iter_iopc_calcp = 3;
        }
        else
        {
            iter_iopc_calcp = 5;
        }

        if (iter_gs)
        {
            iter_gs = !iter_halfe && !iter_ci;
        }

        // determine the self-consistency criterion

        // iter_scfcrt is machine-precision dependent

        iter_scfcrt = 1e-4;

        // increase precision for evrything except normal ground-state calculations

        if (strstr(molkst_c_mp_keywrd_, "PRECISE") || molkst_c_mp_nopen_ != molkst_c_mp_nclose_)
        {
            iter_scfcrt = iter_scfcrt * 0.01;
        }
        iter_scfcrt = max(iter_scfcrt, 1.0e-12);

        if (strstr(molkst_c_mp_keywrd_, "SCFCRT"))
        {
            // something with iter_scfcrt
        }

        if (strstr(molkst_c_mp_keywrd_, "RELSCF"))
        {
            // something with iter_scfcrt
        }

        if (molkst_c_mp_id_ == 3)
        {
            iter_scfcrt = iter_scfcrt * molkst_c_mp_numat_ / 20.0;
        }
        if (iter_scfcrt < 1e-12)
        {
            printf("THERE IS A RISK OF INFINITE LOOPING WITH THE SCFCRT LESS THAN 1.D-12\n");
        }
    }
    else if ((molkst_c_mp_nscf_ > 0) && !molkst_c_mp_uhf_)
    {
        int one = 1;
        dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_pa_, &one, common_arrays_c_mp_pb_, &one);
        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            common_arrays_c_mp_p_[i] = 2.0 * common_arrays_c_mp_pa_[i];
        }
    }

    makea = 1;
    makeb = 1;
    iemin = 0;
    iemax = 0;
    if (iter_irrr != 5)
    {
        if (molkst_c_mp_uhf_)
        {
            int one = 1;
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_pa_, &one, iter_c_mp_pold_, &one);
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_pb_, &one, iter_c_mp_pbold_, &one);
            for (int i = 0; i < molkst_c_mp_norbs_; i++)
            {
                int j = ((i + 1) * (i + 2)) / 2 - 1;
                iter_c_mp_pold2_[i] = common_arrays_c_mp_pa_[j];
                iter_c_mp_pbold2_[i] = common_arrays_c_mp_pb_[j];
            }
        }
        else
        {
            int one = 1;
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_p_, &one, iter_c_mp_pold_, &one);
            for (int i = 0; i < molkst_c_mp_norbs_; i++)
            {
                int j = ((i + 1) * (i + 2)) / 2 - 1;
                iter_c_mp_pold2_[i] = common_arrays_c_mp_p_[j];
            }
        }
    }
    iter_camkin = strstr(molkst_c_mp_keywrd_, "KING") || strstr(molkst_c_mp_keywrd_, "CAMP");

    // turn off iter_shift if not a full scf

    if (!fulscf)
    {
        iter_shift = 0.0;
    }
    if (iter_newdg)
    {
        iter_newdg = fabs(iter_bshift) < 0.001;
    }
    if (molkst_c_mp_last_ == 1)
    {
        iter_newdg = 0;
    }

    // SELF-CONSISTENCY CRITERIA: SELCON IS IN KCAL/MOL, iter_pltest IS
    // A LESS IMPORTANT TEST TO MAKE SURE THAT THE SELCON TEST IS NOT
    // PASSED 'BY ACCIDENT'
    //                            IF GNORM IS LARGE, MAKE SELCON BIGGER

    selcon = iter_scfcrt;

    // LET SELCON BE DETERMINED BY iter_scfcrt AND GNORM, BUT IN NO CASE
    // CAN IT BE MORE THAN 100*SELCON OR 0.1

    if (iter_gs)
    {
        selcon = min(min(iter_scfcrt * 100.0, 0.1), max(iter_scfcrt * pow(molkst_c_mp_gnorm_, 3) * pow(10, -molkst_c_mp_id_ * 3), iter_scfcrt));
    }
    iter_pltest = 0.05 * sqrt(fabs(selcon));

    // SOMETIMES HEAT GOES SCF BUT DENSITY IS STILL FLUCTUATING IN UHF
    // IN WHICH CASE PAY LESS ATTENTION TO DENSITY MATRIX

    if ((molkst_c_mp_nalpha_ != molkst_c_mp_nbeta_) && molkst_c_mp_uhf_)
    {
        iter_pltest = 0.001;
    }

    if (iter_prt1el)
    {
        printf("ONE-ELECTRON MATRIX AT ENTRANCE TO ITER \n");
        vecprt_(common_arrays_c_mp_h_, &molkst_c_mp_norbs_);
    }

    iredy = 1;
l_180:
    // 180 here
    niter = 0;
    frst = 1;
    if (iter_camkin)
    {
        modea = 1;
        modeb = 1;
    }
    else
    {
        modea = 0;
        modeb = 0;
    }

    bfrst = 1;

    // START THE SCF LOOP HERE
l_250:
    // 250 here
    iter_incitr = (modea != 3) && (modeb != 3);
    if (iter_incitr)
    {
        niter += 1;
    }
    if (iter_timitr)
    {
        titer = 0.0;
        if (niter > 1)
        {
            printf("time for iteration blabla wall clock seconds\n");
            iter_titer0 = titer;
        }
    }
    if ((niter > iter_itrmax - 10) && !iter_allcon)
    {
        // SWITCH ON ALL CONVERGERS

        iter_okpuly = 1;
        iter_camkin = !iter_halfe;
        if (iter_itrmax > 2)
        {
            printf("ALL CONVERGERS ARE NOW FORCED ON \n"
                   "iter_shift=10, PULAY ON, CAMP-KING ON AND ITERATION COUNTER RESET\n");
        }
        iter_allcon = 1;
        iter_bshift = 4.44;
        iredy = -4;
        eold = 100.0;
        iter_newdg = 0;
        goto l_180;
    }

    // MAKE THE ALPHA FOCK MATRIX

    if ((fabs(iter_shift) > 1.0e-10) && (iter_bshift != 0.0))
    {
        l = 0;
        if (niter > 1)
        {
            if (iter_newdg && !(iter_halfe || iter_camkin))
            {
                // iter_shift WILL APPLY TO THE VIRTUAL ENERGY LEVELS USED IN THE
                // PSEUDODIAGONALIIZATION. IF DIFF IS -VE, GOOD, THEN LOWER THE
                // HOMO-LUMO GAP BY 0.1EV, OTHERWISE INCREASE IT.
                if (diff > 0.0)
                {
                    iter_shift = 1.0;
                    // IF THE PSEUDODIAGONALIZATION APPROXIMATION -- THAT THE WAVEFUNCTION
                    // IS ALMOST STABLE -- IS INVALID, TURN OFF iter_newdg
                    if (diff > 1)
                    {
                        iter_newdg = 0;
                    }
                }
                else
                {
                    iter_shift = -0.1;
                }
            }
            else
            {
                iter_shift = iter_ten + common_arrays_c_mp_eigs_[ihomo] - common_arrays_c_mp_eigs_[ihomo - 1] + iter_shift;
            }

            if (diff > 0.0)
            {
                if (iter_shift > 4.0)
                {
                    iter_shfmax = 4.5;
                }
                if (iter_shift > iter_shfmax)
                {
                    iter_shfmax = max(iter_shfmax - 0.5, 0);
                }
            }

            // IF SYSTEM GOES UNSTABLE, LIMIT iter_shift TO THE RANGE -INFINITY - iter_shfmax
            // BUT IF SYSTEM IS STABLE, LIMIT iter_shift TO THE RANGE -INFINITY - +20

            iter_shift = max(-20.0, min(iter_shfmax, iter_shift));
            if (fabs(iter_shift < iter_shfmax) < 1e-5)
            {
                iter_shfmax += 0.01;
            }

            // THE CAMP-KING AND PULAY CONVERGES NEED A CONSTANT iter_shift.
            // IF THE iter_shift IS ALLOWED TO VARY, THESE CONVERGERS WILL NOT
            // WORK PROPERLY.

            if (iter_okpuly || (fabs(iter_bshift - 4.44) < 1e-5))
            {
                iter_shift = -8.0;
                if (iter_newdg)
                {
                    iter_shift = 0.0;
                }
            }

            if (molkst_c_mp_uhf_)
            {
                if (iter_newdg && !(iter_halfe || iter_camkin))
                {
                    iter_shiftb = iter_ten - iter_tenold;
                }
                else
                {
                    iter_shiftb = iter_ten + common_arrays_c_mp_eigs_[ihomob] - common_arrays_c_mp_eigs_[ihomob - 1] + iter_shiftb;
                }
                if (diff > 0.0)
                {
                    iter_shiftb = min(4.0, iter_shiftb);
                }
                iter_shiftb = max(-20.0, min(iter_shfmax, iter_shiftb));
                if (iter_okpuly || (fabs(iter_bshift - 4.44) < 1e-5))
                {
                    iter_shiftb = -8.0;
                    if (iter_newdg)
                    {
                        iter_shiftb = 0.0;
                    }
                }
                for (int i = ihomob; i < molkst_c_mp_norbs_; i++)
                {
                    common_arrays_c_mp_eigb_[i] += iter_shiftb;
                }
            }
        }

        iter_tenold = iter_ten;

        if (iter_pl > iter_plchek)
        {
            iter_shftbo = iter_shiftb;
            iter_shfto = iter_shift;
        }
        else
        {
            iter_shiftb = iter_shftbo;
            iter_shift = iter_shfto;
        }

        if (molkst_c_mp_id_ != 0)
        {
            iter_shift = 0.0;
        }

        for (int i = ihomo; i < molkst_c_mp_norbs_; i++)
        {
            common_arrays_c_mp_eigs_[i] += iter_shift;
        }

        if (molkst_c_mp_id_ != 0)
        {
            iter_shift = -80.0;
        }

        if (molkst_c_mp_lxfac_)
        {
            iter_shift = 0.0;
        }

        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            common_arrays_c_mp_f_[i] = common_arrays_c_mp_h_[i] + iter_shift * common_arrays_c_mp_pa_[i] + 1e-16 * (i + 1);
        }

        for (int i = 0; i < molkst_c_mp_norbs_; i++)
        {
            int j = ((i + 1) * (i + 2)) / 2 - 1;
            common_arrays_c_mp_f_[j] -= iter_shift;
        }
    }
    else if ((molkst_c_mp_last_ == 0) && (niter < 2) && fulscf)
    {
        iter_random = 0.001;
        glow = glow || (molkst_c_mp_gnorm_ < 2.0) && (molkst_c_mp_gnorm_ > 1e-9);
        if (glow)
        {
            iter_random = 0.0;
        }
        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            iter_random *= -1; // GBR: This sounds strange. Could iter_random variable be placed out of the loop?
            common_arrays_c_mp_f_[i] = common_arrays_c_mp_h_[i] + iter_random;
        }
    }
    else
    {
        int one = 1;
        dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_h_, &one, common_arrays_c_mp_f_, &one);
    }
l_320:
    // 320 continue
    if (molkst_c_mp_id_ != 0)
    {
        int two = 2;
        fock2(common_arrays_c_mp_f_, common_arrays_c_mp_p_, common_arrays_c_mp_pa_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_wk_, molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, 2, 
              ifact, i1fact, ptot2, &icalcn_fock2, &ione_fock2, jindex_fock2, &lid_fock2);
        //fock2_(common_arrays_c_mp_f_, common_arrays_c_mp_p_, common_arrays_c_mp_c_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_wk_, &molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, &two);
    }
    else
    {
        int two = 2;
        if (dd == 0)
        {
            dd = 1;
        }
        fock2(common_arrays_c_mp_f_, common_arrays_c_mp_p_, common_arrays_c_mp_pa_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, 2,
              ifact, i1fact, ptot2, &icalcn_fock2, &ione_fock2, jindex_fock2, &lid_fock2);
        //fock2_(common_arrays_c_mp_f_, common_arrays_c_mp_p_, common_arrays_c_mp_pa_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, &molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, &two);
    }

    if (molkst_c_mp_uhf_)
    {
        if (iter_shiftb != 0.0)
        {
            l = 0;
            for (int i = 0; i < molkst_c_mp_norbs_; i++)
            {
                for (int j = l; j < i + l; j++)
                {
                    common_arrays_c_mp_fb_[j] = common_arrays_c_mp_h_[j] + iter_shiftb * common_arrays_c_mp_pb_[j];
                    l += i;
                }
                common_arrays_c_mp_fb_[l] -= iter_shiftb;
            }
        }
        else if (/*rand */ (molkst_c_mp_last_ == 0) && (niter < 2) && fulscf)
        {
            iter_random = 0.001;
            if (glow)
            {
                iter_random = 0.0;
            }
            for (int i = 0; i < molkst_c_mp_mpack_; i++)
            {
                iter_random *= -1;
                common_arrays_c_mp_fb_[i] = common_arrays_c_mp_h_[i] + iter_random;
            }
        }
        else
        {
            int one = 1;
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_h_, &one, common_arrays_c_mp_fb_, &one);
        }
        if (molkst_c_mp_id_ != 0)
        {
            int two = 2;
            fock2(common_arrays_c_mp_fb_, common_arrays_c_mp_p_, common_arrays_c_mp_pb_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_wk_, molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, 2, 
                  ifact, i1fact, ptot2, &icalcn_fock2, &ione_fock2, jindex_fock2, &lid_fock2);
            //fock2_(common_arrays_c_mp_fb_, common_arrays_c_mp_p_, common_arrays_c_mp_pb_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_wk_, &molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, &two);
        }
        else
        {
            int two = 2;
            fock2(common_arrays_c_mp_fb_, common_arrays_c_mp_p_, common_arrays_c_mp_pb_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, 2, 
                  ifact, i1fact, ptot2, &icalcn_fock2, &ione_fock2, jindex_fock2, &lid_fock2);
            //fock2_(common_arrays_c_mp_fb_, common_arrays_c_mp_p_, common_arrays_c_mp_pb_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, common_arrays_c_mp_w_, &molkst_c_mp_numat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, &two);
        }
    }
    if (!fulscf)
        goto l_600;
    // CODE THE FOLLOWING LINE IN PROPERLY SOMETIME
    // THIS OPERATION IS BELIEVED TO GIVE RISE TO A BETTER FOCK MATRIX
    // THAN THE CONVENTIONAL GUESS.

    if (iter_irrr == 0)
    {
        for (int i = 0; i < molkst_c_mp_norbs_; i++)
        {
            int j = ((i + 1) * (i + 2)) / 2 - 1;
            common_arrays_c_mp_f_[j] *= 0.5;
        }
    }
    iter_irrr = 2;
    // CALCULATE THE ENERGY IN KCAL/MOLE

    if (niter >= iter_itrmax)
    {
        if ((diff < 1e-3) && (iter_pl < 1e-4) && (!iter_force))
        {
            if (fabs(iter_shift) < 1e-10)
            {
                printf("UNABLE TO ACHIEVE SELF-CONSISTENCE JOB CONTINUING \n");
            }
            iter_incitr = 1;
            getout = 1;
            goto l_410;
        }
        if (iter_minprt)
        {
            printf("UNABLE TO ACHIEVE SELF-CONSISTENCE \n");
        }
        printf("DELTAE= %lf DELTAP= %lf \n", diff, iter_pl);
        molkst_c_mp_iflepo_ = 9;
        molkst_c_mp_iscf_ = 2;
        char *txt = "UNABLE TO ACHIEVE SELF-CONSISTENCE";
        mopend_(txt);
    }
    *ee = helect(molkst_c_mp_norbs_, common_arrays_c_mp_pa_, common_arrays_c_mp_h_, common_arrays_c_mp_f_);
    //printf("ee: %lf\n", *ee);
    if (molkst_c_mp_uhf_)
    {
        *ee = *ee + helect(molkst_c_mp_norbs_, common_arrays_c_mp_pb_, common_arrays_c_mp_h_, common_arrays_c_mp_fb_);
    }
    else
    {
        *ee = *ee * 2.0;
    }
    if (iter_capps)
    {
        *ee = *ee + capcor_(common_arrays_c_mp_nat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, common_arrays_c_mp_p_, common_arrays_c_mp_h_);
    }
    if (molkst_c_mp_uhf_)
    {
        if (iter_bshift != 0.0)
        {
            if (molkst_c_mp_nalpha_open_ > molkst_c_mp_nalpha_)
            {
                iter_scorr = iter_shift * (molkst_c_mp_nalpha_open_ - molkst_c_mp_nalpha_) * iter_enrgy * 0.5 * (molkst_c_mp_fract_ * (1.0 - molkst_c_mp_fract_));
            }
            else
            {
                iter_scorr = iter_shift * (molkst_c_mp_nbeta_open_ - molkst_c_mp_nbeta_) * iter_enrgy * 0.5 * (molkst_c_mp_fract_ * (1.0 - molkst_c_mp_fract_));
            }
        }
    }
    else
    {
        if (iter_bshift != 0.0)
        {
            iter_scorr = iter_shift * (molkst_c_mp_nopen_ - molkst_c_mp_nclose_) * iter_enrgy * 0.25 * (molkst_c_mp_fract_ * (2.0 - molkst_c_mp_fract_));
        }
    }
    escf = (*ee + molkst_c_mp_enuclr_) * iter_enrgy + molkst_c_mp_atheat_ + iter_scorr;
    getout = 0;
l_410:
    if (iter_incitr)
    {
        if (getout)
            goto l_470;
        diff = escf - eold;
        if (diff > 0)
        {
            iter_ten -= 1.0;
        }
        else
        {
            iter_ten = iter_ten * 0.975 + 0.05;
        }
        sellim = max(selcon, 1e-15 * max(fabs(*ee), 1.0));
        // SCF TEST: CHANGE IN HEAT OF FORMATION IN KCAL/MOL SHOULD BE
        // LESS THAN SELLIM.  THE OTHER TESTS ARE SAFETY MEASURES
        if (!(niter > 4 && ((iter_pl == 0.0) || (iter_pl < iter_pltest) && (fabs(diff) < sellim)) && ready))
        {
            goto l_490;
        }
        // SELF-CONSISTENCY TEST, EXIT MODE FROM ITERATIONS
    l_470:
        if (fabs(iter_shift) < 1e-10)
        {
            goto l_600;
        }
        iter_shift = 0.0;
        iter_shiftb = 0.0;
        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            common_arrays_c_mp_f_[i] = common_arrays_c_mp_h_[i];
        }
        makea = 1;
        makeb = 1;
        goto l_320;
    l_490:
        if (molkst_c_mp_limscf_ && (molkst_c_mp_emin_ != 0.0) && !(iter_ci || iter_halfe))
        {
            // THE FOLLOWING TESTS ARE INTENDED TO ALLOW A FAST EXIT FROM ITER
            // IF THE RESULT IS 'GOOD ENOUGH' FOR THE CURRENT STEP IN THE GEOMETRY
            // OPTIMIZATION
            if (escf < molkst_c_mp_emin_)
            {
                // THE ENERGY IS LOWER THAN THE PREVIOUS MINIMUM.  NOW CHECK THAT
                // IT IS CONSISTENTLY LOWER.

                iemax = 0;
                iemin = min(5, iemin + 1);
                for (int i = 0; i < iemin - 1; i++)
                {
                    escf0[i] = escf0[i + 1];
                }
                escf0[iemin - 1] = escf;

                // IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 5%
                // OF THE ENERGY GAIN FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
                // MINIMUM.

                if (iemin > 3)
                {
                    for (int i = 1; i < iemin; i++)
                    {
                        if (fabs(escf0[i] - escf0[i - 1]) > 0.05 * (molkst_c_mp_emin_ - escf))
                        {
                            goto l_540;
                        }
                    }
                    // IS GOOD ENOUGH -- RAPID EXIT
                    iter_incitr = 1;
                    getout = 1;
                    goto l_410;
                }
            }
            else
            {
                // THE ENERGY HAS RISEN ABOVE THAT OF THE PREVIOUS MINIMUM.
                // WE NEED TO CHECK WHETHER THIS IS A FLUKE OR IS THIS REALLY
                // A BAD GEOMETRY.
                iemin = 0;
                iemax = min(5, iemax + 1);
                for (int i = 0; i < iemax - 1; i++)
                {
                    escf0[i] = escf0[i + 1];
                }
                escf0[iemax - 1] = escf;
                // IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 5%
                // OF THE ENERGY LOST FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
                // MINIMUM.

                if (iemax > 3)
                {
                    for (int i = 1; i < iemax; i++)
                    {
                        if (fabs(escf0[i] - escf0[i - 1]) > 0.05 * (escf - molkst_c_mp_emin_))
                        {
                            goto l_540;
                        }
                        iter_incitr = 1;
                        getout = 1;
                        goto l_410;
                    }
                }
            }
        }
    l_540:
        ready = (iredy > 0) && ((fabs(diff) < sellim * 10.0) || (iter_pl == 0.0));
        iredy += 1;
    }
    printf("ITERATION %d PLS=%lf ENERGY=%lf, DELTAE=%lf \n", niter, iter_pl, escf, diff);
    if (iter_incitr)
    {
        eold = escf;
    }
    // INVOKE THE CAMP-KING CONVERGER

    if ((niter > 2) && (iter_camkin) && makea)
    {
        if (iter_c_mp_vec_ai_ == NULL)
        {
            iter_c_mp_vec_ai_ = (double *)malloc(molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double));
            iter_c_mp_fock_ai_ = (double *)malloc(molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double));
            iter_c_mp_p_ai_ = (double *)malloc(molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double));
            iter_c_mp_h_ai_ = (double *)malloc(molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double));
            iter_c_mp_vecl_ai_ = (double *)malloc(molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double));
            // TODO: make it parallel
            for (int i = 0; i < molkst_c_mp_norbs_ * molkst_c_mp_norbs_; i++)
            {
                iter_c_mp_vec_ai_[i] = 0.0;
                iter_c_mp_fock_ai_[i] = 0.0;
                iter_c_mp_p_ai_[i] = 0.0;
                iter_c_mp_h_ai_[i] = 0.0;
                iter_c_mp_vecl_ai_[i] = 0.0;
            }
        }
        int norbs_na1el = molkst_c_mp_norbs_ - iter_na1el;
        double escf_enrgy = escf / iter_enrgy;
        interp_(&iter_na1el, &norbs_na1el, &modea, &escf_enrgy, common_arrays_c_mp_f_, common_arrays_c_mp_c_, theta, iter_c_mp_vec_ai_, iter_c_mp_fock_ai_, iter_c_mp_p_ai_, iter_c_mp_h_ai_, iter_c_mp_vecl_ai_, &eold_alpha);
    }
    makeb = 0;
    if (modea != 3)
    {
        makeb = 1;
        if (iter_newdg)
        {
            // INVOKE PULAY'S CONVERGER
            if (iter_okpuly && makea && (iredy > 1))
            {
                int mpack_6 = 6 * molkst_c_mp_mpack_;
                pulay_(common_arrays_c_mp_f_, common_arrays_c_mp_pa_, &molkst_c_mp_norbs_, iter_c_mp_pold_, iter_c_mp_pold2_, iter_c_mp_pold3_, &iter_jalp, &iter_ialp, &mpack_6, &frst, &iter_pl);
            }

            // DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT
            // WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S

            if (iter_halfe || iter_camkin)
            {
                eigenvectors_lapack_(common_arrays_c_mp_c_, common_arrays_c_mp_f_, common_arrays_c_mp_eigs_, &molkst_c_mp_norbs_);
            }
            else
            {
                // if cpu
                diag_for_gpu_(common_arrays_c_mp_f_, common_arrays_c_mp_c_, &iter_na1el, common_arrays_c_mp_eigs_, &molkst_c_mp_norbs_, &molkst_c_mp_mpack_);
            }
        }
        else
        {
            eigenvectors_lapack_(common_arrays_c_mp_c_, common_arrays_c_mp_f_, common_arrays_c_mp_eigs_, &molkst_c_mp_norbs_);
        }

        if (iter_ifill != 0)
        {
            swap_(common_arrays_c_mp_c_, &molkst_c_mp_norbs_, &molkst_c_mp_norbs_, &iter_na2el, &iter_ifill);
        }

        // CALCULATE THE ALPHA OR RHF DENSITY MATRIX

        if (molkst_c_mp_uhf_)
        {
            double one = 1.0;
            int one_int = 1;
            density_for_gpu_(common_arrays_c_mp_c_, &molkst_c_mp_fract_, &molkst_c_mp_nalpha_, &molkst_c_mp_nalpha_open_, &one, &molkst_c_mp_mpack_, &molkst_c_mp_norbs_, &one_int, common_arrays_c_mp_pa_, &iter_iopc_calcp);
            if ((modea != 3) && !(iter_newdg && iter_okpuly))
            {
                int ii = niter;
                if (iter_camkin)
                {
                    ii = 7;
                }
                cnvg_(common_arrays_c_mp_pa_, iter_c_mp_pold_, iter_c_mp_pold2_, &ii, &iter_pl);
            }
        }
        else
        {
            if (iter_halfe)
            {
                double two = 2.0;
                int one = 1;
                densit_(common_arrays_c_mp_c_, &molkst_c_mp_norbs_, &molkst_c_mp_norbs_, &iter_na2el, &two, &iter_na1el, &molkst_c_mp_fract_, common_arrays_c_mp_p_, &one);
            }
            else
            {
                double two = 2.0;
                int one = 1;
                density_for_gpu_(common_arrays_c_mp_c_, &molkst_c_mp_fract_, &iter_na2el, &iter_na1el, &two, &molkst_c_mp_mpack_, &molkst_c_mp_norbs_, &one, common_arrays_c_mp_p_, &iter_iopc_calcp);
            }

            if ((modea != 3) && !(iter_newdg && iter_okpuly))
            {
                cnvg_(common_arrays_c_mp_p_, iter_c_mp_pold_, iter_c_mp_pold2_, &niter, &iter_pl);
            }
        }
    }

    // UHF-SPECIFIC CODE

    if (molkst_c_mp_uhf_)
    {
        // INVOKE THE CAMP-KING CONVERGER

        if ((niter > 2) && (iter_camkin) && (makeb))
        {
            if (iter_c_mp_vec_bi_ == NULL)
            {
                double all_size = molkst_c_mp_norbs_ * molkst_c_mp_norbs_ * sizeof(double);
                iter_c_mp_vec_bi_ = (double *)malloc(all_size);
                iter_c_mp_fock_bi_ = (double *)malloc(all_size);
                iter_c_mp_p_bi_ = (double *)malloc(all_size);
                iter_c_mp_h_bi_ = (double *)malloc(all_size);
                iter_c_mp_vecl_bi_ = (double *)malloc(all_size);
            }

            int norbs_nb1el = molkst_c_mp_norbs_ - iter_nb1el;
            double escf_enrgy = escf / iter_enrgy;
            interp_(&iter_nb1el, &norbs_nb1el, &modeb, &escf_enrgy, common_arrays_c_mp_fb_, common_arrays_c_mp_cb_, theta, iter_c_mp_vec_bi_, iter_c_mp_fock_bi_, iter_c_mp_p_bi_, iter_c_mp_h_bi_, iter_c_mp_vecl_bi_, &eold_beta);
        }
        makea = 0;
        if (modeb != 3)
        {
            makea = 1;
            if (iter_newdg)
            {
                // INVOKE PULAY'S CONVERGER
                if (iter_okpuly && makeb && (iredy > 1))
                {
                    int mpack_6 = 6 * molkst_c_mp_mpack_;
                    pulay_(common_arrays_c_mp_fb_, common_arrays_c_mp_pb_, &molkst_c_mp_norbs_, iter_c_mp_pbold_, iter_c_mp_pbold2_, iter_c_mp_pbold3_, &iter_jbet, &iter_ibet, &mpack_6, &bfrst, &iter_plb);
                }

                // DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT
                // WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S

                if (iter_halfe || iter_camkin)
                {
                    eigenvectors_lapack_(common_arrays_c_mp_cb_, common_arrays_c_mp_fb_, common_arrays_c_mp_eigb_, &molkst_c_mp_norbs_);
                }
                else
                {
                    // if cpu
                    diag_for_gpu_(common_arrays_c_mp_fb_, common_arrays_c_mp_cb_, &iter_nb1el, common_arrays_c_mp_eigb_, &molkst_c_mp_norbs_, &molkst_c_mp_mpack_);
                }
            }
            else
            {
                eigenvectors_lapack_(common_arrays_c_mp_cb_, common_arrays_c_mp_fb_, common_arrays_c_mp_eigb_, &molkst_c_mp_norbs_);
            }
        }

        // CALCULATE THE BETA DENSITY MATRIX

        double gpu_one = 1.0;
        int int_one = 1;
        density_for_gpu_(common_arrays_c_mp_cb_, &molkst_c_mp_fract_, &molkst_c_mp_nbeta_, &molkst_c_mp_nbeta_open_, &gpu_one, &molkst_c_mp_mpack_, &molkst_c_mp_norbs_, &int_one, common_arrays_c_mp_pb_, &iter_iopc_calcp);
        if (!(iter_newdg && iter_okpuly))
        {
            int ii = niter;
            if (iter_camkin)
            {
                ii = 7;
            }
            cnvg_(common_arrays_c_mp_pb_, iter_c_mp_pbold_, iter_c_mp_pbold2_, &ii, &iter_plb);
        }
    }

    // calculate the total density matrix

    if (molkst_c_mp_uhf_)
    {
        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            common_arrays_c_mp_p_[i] = common_arrays_c_mp_pa_[i] + common_arrays_c_mp_pb_[i];
        }
    }
    else
    {
        for (int i = 0; i < molkst_c_mp_mpack_; i++)
        {
            common_arrays_c_mp_pa_[i] = common_arrays_c_mp_p_[i] * 0.5;
            common_arrays_c_mp_pb_[i] = common_arrays_c_mp_pa_[i];
        }
    }

    if (iter_itrmax < 3)
    {
        return;
    }

    iter_oknewd = (iter_pl < sellim) || iter_oknewd;
    iter_newdg = (iter_pl < iter_trans) && iter_oknewd || iter_newdg;
    if (iter_pl < iter_trans * 0.3333)
    {
        iter_oknewd = 1;
    }
    goto l_250;

    // END THE SCF LOOP HERE
    // NOW CALCULATE THE ELECTRONIC ENERGY

    // SELF-CONSISTENCE ACHIEVED.

l_600:
    *ee = helect(molkst_c_mp_norbs_, common_arrays_c_mp_pa_, common_arrays_c_mp_h_, common_arrays_c_mp_f_);
    if (molkst_c_mp_uhf_)
    {
        *ee = *ee + helect(molkst_c_mp_norbs_, common_arrays_c_mp_pb_, common_arrays_c_mp_h_, common_arrays_c_mp_fb_);
    }
    else
    {
        *ee = *ee * 2.0;
    }
    if (iter_capps)
    {
        *ee = *ee + capcor_(common_arrays_c_mp_nat_, common_arrays_c_mp_nfirst_, common_arrays_c_mp_nlast_, common_arrays_c_mp_p_, common_arrays_c_mp_h_);
    }
    if ((molkst_c_mp_nscf_ == 0) || (molkst_c_mp_last_ == 1) || iter_ci || iter_halfe)
    {
        // PUT F AND FB INTO POLD IN ORDER TO NOT DESTROY F AND FB
        // AND DO EXACT DIAGONALISATIONS
        int one = 1;
        dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_f_, &one, iter_c_mp_pold_, &one);
        eigenvectors_lapack_(common_arrays_c_mp_c_, iter_c_mp_pold_, common_arrays_c_mp_eigs_, &molkst_c_mp_norbs_);
        if (molkst_c_mp_last_ == 1)
        {
            phase_lock_(common_arrays_c_mp_c_, &molkst_c_mp_norbs_);
        }
        if (molkst_c_mp_uhf_)
        {
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_fb_, &one, iter_c_mp_pold_, &one);
            eigenvectors_lapack_(common_arrays_c_mp_cb_, iter_c_mp_pold_, common_arrays_c_mp_eigb_, &molkst_c_mp_norbs_);
            if (molkst_c_mp_last_ == 1)
            {
                phase_lock_(common_arrays_c_mp_cb_, &molkst_c_mp_norbs_);
            }
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_pa_, &one, iter_c_mp_pold_, &one);
        }
        else
        {
            dcopy_(&molkst_c_mp_mpack_, common_arrays_c_mp_p_, &one, iter_c_mp_pold_, &one);
        }
        if (iter_ci || iter_halfe)
        {
            sum = meci_();
            if (molkst_c_mp_moperr_)
            {
                return;
            }
            *ee = *ee + sum;
            if (iter_prtpl)
            {
                escf = (*ee + molkst_c_mp_enuclr_) * iter_enrgy + molkst_c_mp_atheat_;
            }
        }
    }
    molkst_c_mp_nscf_ += 1;
    if (iter_allcon && (fabs(iter_bshift - 4.44) < 1e-7))
    {
        iter_camkin = 0;
        iter_allcon = 0;
        iter_newdg = 0;
        iter_bshift = -10.0;
        iter_okpuly = 0;
    }

    iter_shift = 1.0;
    escf = (*ee + molkst_c_mp_enuclr_) * iter_enrgy + molkst_c_mp_atheat_;
    if (molkst_c_mp_emin_ == 0.0)
    {
        molkst_c_mp_emin_ = escf;
    }
    else
    {
        molkst_c_mp_emin_ = min(molkst_c_mp_emin_, escf);
    }

    for (int i = 0; i < 81; i++)
    {
        free(ptot2[i]);
    }

    free(ptot2);
    free(ifact);
    free(i1fact);

    return;
}
