#include <math.h>

#include <stdlib.h>

#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

extern int molkst_c_mp_numcal_;
extern int molkst_c_mp_norbs_;
extern int molkst_c_mp_mpack_;
extern int molkst_c_mp_n2elec_;
extern int molkst_c_mp_id_;
extern int molkst_c_mp_numat_;

extern int cosmo_c_mp_useps_;

extern void fockdorbs_(int *ia, int *ib, int *ja, int *jb, double *f, double *p, double *ptot, double *w, int *kr, int *ifact);
extern void fockdorbs(int ia, int ib, int ja, int jb, double *f, double *p, double *ptot, double *w, int *kr, int *ifact);
extern void jab_(int *IA, int *JA, double *PJA, double *PJB, double *W, double *F);
extern void kab_(int *IA, int *JA, double *PK, double *W, double *F);

extern void fock1dorbs_(double *f, double *ptot, double *pa, int *mpack, double *w, int *kr, int *ia, int *ib, int *ilim);


extern void addfck_(double *f, double *ptot);

void fock2(double *f, double *ptot, double *p, double *w, double *wj, double *wk,
           int numat, int *nfirst, int *nlast, int mode, int *ifact, int *i1fact,
           double **ptot2, int *icalcn, int *ione, int *jindex, int *lid)
{

    if (numat == 0)
    {
        if (ptot2)
        {
            for (int i = 0; i < 81; i++)
            {
                free(ptot2[i]);
            }
            free(ptot2);
        }
        if (ifact)
            free(ifact);
        if (i1fact)
            free(i1fact);
        return;
    }

    double pk[16];
    double pja[16];
    double pjb[16];

    double sumdia, sumoff, sum, aa, bb, aj, ak, a;

    int deriv;

    int ij, ji, ik, kl, lk, jl, il, iminus, jk;

    deriv = (numat < 0);
    numat = abs(numat);

    if (*icalcn != molkst_c_mp_numcal_)
    {
        *icalcn = molkst_c_mp_numcal_;

        // SET UP ARRAY OF LOWER HALF TRIANGLE INDICES (PASCAL'S TRIANGLE)

        for (int i = 0; i < molkst_c_mp_norbs_; i++)
        {
            ifact[i] = (i * (i + 1)) / 2;
            i1fact[i] = ifact[i] + i + 1;
        }

        // SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
        // INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
        // INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
        //              KINDEX ARE THE INDICES OF THE K-INTEGRALS

        int m = -1;
        for (int i = 1; i <= 4; i++)
        {
            for (int j = 1; j <= 4; j++)
            {
                ij = min(i, j);
                ji = i + j - ij;
                for (int k = 1; k <= 4; k++)
                {
                    ik = min(i, k);
                    for (int l = 1; l <= 4; l++)
                    {
                        m += 1;
                        kl = min(k, l);
                        lk = k + l - kl;
                        jl = min(j, l);
                        jindex[m] = (ifact[ji - 1] + ij) * 10 + ifact[lk - 1] + kl - 10;
                    }
                }
            }
        }
        *lid = (molkst_c_mp_id_ == 0);
        *ione = *lid;

        // END OF INITIALIZATION
    }

    // START OF MNDO, AM1 OR PM3 OPTION

    for (int i = 0; i < numat; i++)
    {
        int ia = nfirst[i];
        int ib = nlast[i];
        int m = -1;
        for (int j = ia; j <= ib; j++)
        {
            for (int k = ia; k <= ib; k++)
            {
                m += 1;
                int jk = min(j, k);
                int kj = k + j - jk;
                jk += (kj * (kj - 1)) / 2;
                ptot2[m][i] = ptot[jk - 1];
            }
        }
    }

    int kk = 0;

    for (int ii = 1; ii <= numat; ii++)
    {
        int ia = nfirst[ii - 1];
        int ib = nlast[ii - 1];

        // IF NUMAT=2 THEN WE ARE IN A DERIVATIVE IN A SOLID STATE OR IN A MOLECULE CALCULATION

        if (deriv)
        {
            iminus = ii - 1;
        }
        else
        {
            iminus = ii - *ione;
        }

        for (int jj = 1; jj <= iminus; jj++)
        {
            int ja = nfirst[jj - 1];
            int jb = nlast[jj - 1];

            if (*lid)
            {
                if ((ib - ia >= 6) || (jb - ja >= 6))
                {
                    //fockdorbs_(&ia, &ib, &ja, &jb, f, p, ptot, w, &kk, ifact);
                    fockdorbs(ia, ib, ja , jb, f, p, ptot, w, &kk, ifact);
                }
                else if ((ib - ia >= 3) && (jb - ja) >= 3)
                {
                    // HEAVY-ATOM  - HEAVY-ATOM
                    // EXTRACT COULOMB TERMS

                    for (int i = 0; i < 16; i++)
                    {
                        pja[i] = ptot2[i][ii - 1];
                        pjb[i] = ptot2[i][jj - 1];
                    }

                    jab_(&ia, &ja, pja, pjb, &(w[kk]), f);

                    // EXCHANGE TERMS

                    // EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX

                    // The following loop had been written in a more compact style, but
                    // this had caused problems with the INTEL compiler in RELEASE mode when QuickWin was used

                    int ll = 1;
                    for (int i = ia - 1; i < ib; i++)
                    {
                        int i1 = ifact[i] + ja - 1;
                        for (int j = ll - 1; j < ll + 3; j++)
                        {
                            pk[j] = p[i1];
                            i1 += 1;
                        }
                        ll += 4;
                    }

                    kab_(&ia, &ja, pk, &(w[kk]), f);

                    kk += 100;
                }
                else if ((ib - ia >= 3) && (ja == jb))
                {
                    // LIGHT-ATOM - HEAVY-ATOM

                    // COULOMB TERMS

                    sumdia = 0.0;
                    sumoff = 0.0;

                    int lll = i1fact[ja - 1] - 1;
                    int k = 0;
                    for (int i = 0; i <= 3; i++)
                    {
                        int j1 = ifact[ia + i - 1] + ia - 1;
                        if (i > 0)
                        {
                            for (int j = 1; j <= i; j++)
                            {
                                f[j + j1 - 1] += ptot[lll] * w[j + kk + k - 1];
                                sumoff += ptot[j + j1 - 1] * w[j + kk + k - 1];
                            }
                            k += i;
                            j1 += i;
                        }
                        j1 += 1;
                        k += 1;
                        f[j1 - 1] += ptot[lll] * w[kk + k - 1];
                        sumdia += ptot[j1 - 1] * w[kk + k - 1];
                    }
                    f[lll] += sumoff * 2.0 + sumdia;
////////////////////////////////////////////////////
                    // EXCHANGE TERMS
                    // EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX

                    k = 0;
                    for (int i = ia - 1; i < ib; i++)
                    {
                        int i1 = ifact[i] + ja;
                        sum = 0.0;
                        for (int j = 0; j < ib - ia + 1; j++)
                        {
                            sum += p[ifact[j - 1 + ia] + ja - 1] * w[kk + jindex[j + k] - 1];
                        }
                        k += ib - ia + 1;
                        f[i1 - 1] -= sum;
                    }
                    kk += 10;
                }
                else if ((jb - ja >= 3) && (ia == ib))
                {
                    // HEAVY-ATOM - LIGHT-ATOM

                    // COULOMB TERMS

                    sumdia = 0.0;
                    sumoff = 0.0;
                    int lll = i1fact[ia-1] - 1;
                    int k = 0;
                    for (int i = 0; i <= 3; i++)
                    {
                        int j1 = ifact[ja + i - 1] + ja - 1;
                        if (i > 0)
                        {
                            for (int j = 1; j <= i; j++)
                            {
                                f[j + j1 - 1] += ptot[lll] * w[j + kk + k - 1];
                                sumoff += ptot[j + j1 - 1] * w[j + kk + k - 1];
                            }
                            k += i;
                            j1 += i;
                        }
                        j1 += 1;
                        k += 1;
                        f[j1 - 1] += ptot[lll] * w[kk + k - 1];
                        sumdia += ptot[j1 - 1] * w[kk + k - 1];
                    }
                    f[lll] += sumoff * 2.0 + sumdia;
                    // EXCHANGE TERMS
                    // EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX

                    k = ifact[ia - 1] + ja;
                    int j = 0;
                    for (int i = k - 1; i < k + 3; i++)
                    {
                        sum = 0.0;
                        for (int l = 0; l < 4; l++)
                        {
                            sum += p[l - 1 + k] * w[kk + jindex[l + j] - 1];
                        }
                        j += 4;
                        f[i] -= sum;
                    }
                    kk += 10;
                }
                else if ((jb == ja) && (ia == ib))
                {
                    // LIGHT-ATOM - LIGHT-ATOM

                    int i1 = i1fact[ia - 1];
                    int j1 = i1fact[ja - 1];
                    ij = i1 + ja - ia;
                    f[i1 - 1] += ptot[j1 - 1] * w[kk];
                    f[j1 - 1] += ptot[i1 - 1] * w[kk];
                    f[ij - 1] -= p[ij - 1] * w[kk];
                    kk += 1;
                }
            }
            else
            {
                for (int i = ia; i <= ib; i++)
                {
                    int ka = ifact[i-1];
                    for (int j = ia; j <= i; j++)
                    {
                        int kb = ifact[j-1];
                        ij = ka + j;
                        aa = 2.0;
                        if (i == j)
                        {
                            aa = 1.0;
                        }
                        for (int k = ja; k <= jb; k++)
                        {
                            int kc = ifact[k-1];
                            if (i >= k)
                            {
                                ik = ka + k;
                            }
                            else
                            {
                                ik = 0;
                            }
                            if (j >= k)
                            {
                                jk = kb + k;
                            }
                            else
                            {
                                jk = 0;
                            }

                            for (int l = ja; l <= k; l++)
                            {
                                if (i >= l)
                                {
                                    il = ka + l;
                                }
                                else
                                {
                                    il = 0;
                                }
                                if (j >= l)
                                {
                                    jl = kb + l;
                                }
                                else
                                {
                                    jl = 0;
                                }
                                kl = kc + l;
                                bb = 2.0;
                                if (k == l)
                                {
                                    bb = 1.0;
                                }
                                kk = kk + 1;
                                aj = wj[kk - 1];
                                ak = wk[kk - 1];

                                // A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
                                // ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
                                // AA AND BB ARE CORRECTION FACTORS SINCE
                                // (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
                                // IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIB ORBITALS
                                // I AND J.  SIMILARLY FOR IK ETC.

                                // THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
                                // MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
                                // BUT IT WORKS,

                                if (kl > ij)
                                    continue;
                                if ((i == k) && (aa + bb < 2.1))
                                {
                                    f[ij-1] += aj * ptot[kl-1];
                                }
                                else
                                {
                                    f[ij-1] += bb * aj * ptot[kl-1];
                                    f[kl-1] += aa * aj * ptot[ij-1];
                                    a = ak * aa * bb * 0.25;
                                    if (jl > 0)
                                        f[ik-1] -= a * p[jl-1];
                                    if (jk > 0)
                                        f[il-1] -= a * p[jk-1];
                                    if (jk > 0)
                                        f[jk-1] -= a * p[il-1];
                                    if (jl > 0)
                                        f[jl-1] -= a * p[ik-1];
                                }
                            }
                        }
                    }
                }
            }
        }
        if (mode == 2)
        {
            int i = ((ib - ia + 1) * (ib - ia + 2)) / 2;
            fock1dorbs_(f, ptot, p, &molkst_c_mp_mpack_, &(w[kk]), &kk, &ia, &ib, &i);
        }
    }

    if (cosmo_c_mp_useps_)
    {
        addfck_(f, ptot);
    }
    return;
}