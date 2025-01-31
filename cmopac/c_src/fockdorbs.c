

void fockdorbs(int ia, int ib, int ja, int jb, double *f, double *p, double *ptot, double *w, int *kr, int *ifact){
    if(ia > ja){
        for(int i = ia; i <= ib; i++){
            int ka = ifact[i-1];
            double aa = 2.0;
            for(int j = ia; j <= i; j++){
                if(i == j){
                    aa = 1.0;
                }
                int kb = ifact[j-1];
                int ij = ka + j;
                for(int k = ja; k <= jb; k++){
                    int kc = ifact[k-1];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for(int l = ja; l <= k; l++){
                        if(k == l){
                            bb = 1.0;
                        }
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        *kr += 1;
                        double a = w[*kr - 1];
                        f[ij - 1] += bb * a * ptot[kl - 1]; 
                        f[kl - 1] += aa * a * ptot[ij - 1];
                        a *= aa * bb * 0.25;
                        f[ik - 1] -= a * p[jl - 1];
                        f[il - 1] -= a * p[jk - 1];
                        f[jk - 1] -= a * p[il - 1];
                        f[jl - 1] -= a * p[ik - 1];
                    }
                }
            }
        }
    }
    else {
        int kref = *kr;
        int nn = jb - ja;
        nn = ((nn) * (nn + 1)) / 2;
        int n1 = 0;
        for(int i = ja; i <= jb; i++){
            int ka = ifact[i - 1];
            double aa = 2.0;
            for(int j = ja; j <= i; j++){
                n1 += 1;
                if(i == j){
                    aa = 1.0;
                }
                int kb = ifact[j - 1];
                int ij = ka + j;
                int n2 = 0;
                for(int k = ia; k <= ib; k++){
                    int kc = ifact[k - 1];
                    int ik = ka + k;
                    int jk = kb + k;
                    double bb = 2.0;
                    for(int l = ia; l <= k; l++){
                        n2 += 1;
                        if(k == l){
                            bb = 1.0;
                        }
                        int il = ka + l;
                        int jl = kb + l;
                        int kl = kc + l;
                        *kr += 1;
                        double a = w[kref + (n2-1)*nn + n1];
                        f[ij - 1] += bb * a * ptot[kl - 1];
                        f[kl - 1] += aa * a * ptot[ij - 1];
                        a *= aa * bb * 0.25;
                        f[ik - 1] -= a * p[jl - 1];
                        f[il - 1] -= a * p[jk - 1];
                        f[jk - 1] -= a * p[il - 1];
                        f[jl - 1] -= a * p[ik - 1];
                    }
                }
            }
        }
    }
}