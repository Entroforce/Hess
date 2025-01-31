double helect(int n, double *p, double *h, double *f){
    double ed = 0.0;
    double ee = 0.0;
    int k = -1;
    int nn = n + 1;
    for(int i = 2; i <= nn; i++) {
        k++;
        int jj = i - 1;
        ed += p[k] * (h[k] + f[k]);
        if(i == nn)
            continue;
        if(jj > 0){
            for(int j = k+1; j <= jj+k; j++){
                ee += p[j] * (h[j] + f[j]);
            }
            k += jj;
        }
    }
    ee += 0.5*ed;
    return ee;
}