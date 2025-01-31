#include <stdio.h>
#include <string.h>

double run_mopac(unsigned char *jobnam, int* ck);
int main(int argc, char *argv[]){
    unsigned char jobnam[241];
    for(int i = 0; i < strlen(argv[1]); i++){
        jobnam[i] = argv[1][i];
    }
    for(int i = strlen(argv[1]); i < 240; i++){
        jobnam[i]= ' ';
    }
    jobnam[240] = '\0';
    int one = 1;
    double energy = run_mopac(jobnam, &one);
    printf("end energy: %lf \n", energy);
    return 0;
}