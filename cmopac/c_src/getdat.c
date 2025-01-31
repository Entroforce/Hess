#include <stdio.h>
#include <stdbool.h> 
#include <string.h>
#include "config.h"


extern int molkst_c_mp_natoms_;
extern int molkst_c_mp_ncomments_;

int ReadLine(char *buff, int size, FILE *fp){
    buff[0] = '\0';
    buff[size - 1] = '\0'; /* mark end of buffer */
    char *tmp;

    if(fgets(buff, size, fp) == NULL){
        *buff = '\0';
        return false;
    }
    else {
        /* remove newline */
        if((tmp = strrchr(buff, '\n')) != NULL){
            *tmp = '\0';
        }
    }

    return true;
}

bool getdat(Config *config) {
    int from_data_set = 7;
    int i, j, io_stat, l, nlines, ncomments;
    bool exists, arc_file = false, comments = true;
    char text[90], line1[1000], num1, num2;
    char tmp_comments[120];

    molkst_c_mp_natoms_ = 1;
    // check if an argument is proper filename
    printf("%s\n", config->args);
    char *dot = strrchr(config->args,'.');
    if(dot && !(strcmp(dot, ".mop") && strcmp(dot, ".dat") && strcmp(dot, ".arc"))){
        printf("Filename ends with .mop or .dat or .arc\n");
        if(!strcmp(dot, "arc")){
            arc_file = true;
        }
        config->input = fopen(config->args, "r");

        if(config->input == NULL){
            return false;
        }

        nlines = 0;
        ncomments = 0;
        
        char c;
        ReadLine(config->keywrd, 1000, config->input);
        int count = 2;
        for (c = getc(config->input); c != EOF; c = getc(config->input)){
            if( c == '\n'){
                count += 1;
            }
        }
        config->natoms = count + 200;
        config->maxatoms = config->natoms;
        printf("%s\n", config->keywrd);
        printf("%d\n", count);

    } else {
        return false;
    }

    return true;
}