/**************************************************************************
 * This file is part of the Hess project
 * Copyright (C) 2024 Entroforce LLC
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 **************************************************************************/

#include "hess.h"
#include "hess_se.h"
#include <string.h>
#include <stdlib.h>

typedef struct {
  char** mol;
  char** out;
  char** optimize;
  int n;
} SEParams;


int SECmdParse(int argc, char *argv[], SEParams* params) {
    const char *prot_mol = "-p";
    const char *lig_mol = "-l";
    const char *prot_opt = "-popt";
    const char *lig_opt = "-lopt";
    const char *prot_out = "-pout";
    const char *lig_out = "-lout";


  params->n = 2;
  params->mol = (char **) malloc(params->n * sizeof(char *));
  params->out = (char **) malloc(params->n * sizeof(char *));
  params->optimize = (char **) malloc(params->n * sizeof(char *));

  const char *help_message =
          "Input:\n"
          "-p arg         path to the protein\n"
          "-l arg         path to the ligand\n"
          "Minimization options:\n"
          "-popt arg      parameters for the protein optimization space ('n' for no optimization; 'y' for all molecule optimization; 'h' for only hydrogens optimization)\n"
          "-lopt arg      parameters for the ligand optimization space ('n' for no optimization; 'y' for all molecule optimization; 'h' for only hydrogens optimization)\n"
          "Output: \n"
          "-pout arg      output protein file name, sdf format is used\n"
          "-lout arg      output ligand file name, sdf format is used\n"
          ;
  int i = 0;
  while (i < argc) {
    const char* curr = argv[i];
    if (strcmp(curr, "-h") == 0 || strcmp(curr, "--help") == 0) {
      fprintf(stderr, "%s", help_message);
      return -1;
    }
    
    if ((strlen(curr) && curr[0] != '-') || (strlen(curr) > 1 && curr[0] != '-' && curr[1] != '-')){
        i++;
        continue;
    }
    
    printf("%s\n", curr);
      
    if (strcmp(curr, prot_mol) == 0 && i + 1 < argc){
        params->mol[0] = strdup(argv[i + 1]);
        i++;
    }
    
    else if (strcmp(curr, lig_mol) == 0 && i + 1 < argc) {
        params->mol[1] = strdup(argv[i + 1]);
        i++; 
    }
    
    else if (strcmp(curr, prot_out) == 0 && i + 1 < argc){
        params->out[0] = strdup(argv[i + 1]);
        i++;
    }
    
    else if (strcmp(curr, lig_out) == 0 && i + 1 < argc){
        params->out[1] = strdup(argv[i + 1]);
        i++;
    }

    else if (strcmp(curr, prot_opt) == 0 && i + 1 < argc){
        params->optimize[0] = strdup(argv[i + 1]);
        if(strcmp(params->optimize[0], "y") != 0 && strcmp(params->optimize[0], "n") != 0 && strcmp(params->optimize[0], "h") != 0) {
            fprintf(stderr, "Unrecognized option %s.'\n\nCorrect usage: %s\n\n", params->optimize[0], help_message);
            return -1;
        }
        i++;
    }
    
    else if (strcmp(curr, lig_opt) == 0 && i + 1 < argc){
        params->optimize[1] = strdup(argv[i + 1]);
        if(strcmp(params->optimize[1], "y") != 0 && strcmp(params->optimize[1], "n") != 0 && strcmp(params->optimize[1], "h") != 0) {
            fprintf(stderr, "Unrecognized option %s.'\n\nCorrect usage: %s\n\n", params->optimize[1], help_message);
            return -1;
        }
        i++;
    }
    
    else {
      fprintf(stderr,
              "Unrecognized option %s.\n\nCorrect usage: %s\n\n", curr, help_message);
      return -1;
    }
    
    i++;
  }

  return 0;
}

int main(int argc, char** argv) {

    SEParams params;
    if (SECmdParse(argc, argv, &params) != 0) {
        fprintf(stderr, "%s\n", hessGetLastError());
        return -1;
    }
    
    hessInit();
    hessSetStream(stdout);
    
    void* parser = hessCreateParser();
    
 //   void * se_complex = hessSECreate();
    void * se_prot = hessSECreate();
   // void * se_lig = hessSECreate();
    
    int *indices = (int *)malloc(params.n * sizeof(int));
    int res;
    
    int optModeComplex = 0;
    int optModeProt = 0;
    int optModeLig = 0;
    
    // read protein
    
   if(strcmp(params.optimize[0], "y") == 0) {
        optModeProt = 1;
    } else if(strcmp(params.optimize[0], "h") == 0) {
        optModeProt = 2;
    }
    
    void* atomsProt = hessLoadMolecule(parser, params.mol[0]);
    if (atomsProt == NULL) {
        hessDestroy(parser);
        fprintf(stderr, "Error reading the receptor file. ");
        fprintf(stderr, "%s\n", hessGetLastError());
        return -1;
    }

    hessAssignBonds(atomsProt);
    hessProtonate(atomsProt);
    
    
    int idxProt = 0;    
    res = hessSEAddMol(se_prot, atomsProt, optModeProt, &idxProt);
    
    // read ligand
    
    // if(strcmp(params.optimize[1], "y") == 0) {
    //     optModeLig = 1;
    // } else if(strcmp(params.optimize[1], "h") == 0) {
    //     optModeLig = 2;
    // }
    
    // void* atomsLig = hessLoadMolecule(parser, params.mol[1]);
    // if (atomsLig == NULL) {
    //     hessDestroy(parser);
    //     fprintf(stderr, "Error reading the receptor file. ");
    //     fprintf(stderr, "%s\n", hessGetLastError());
    //     return -1;
    // }

//    hessAssignBonds(atomsLig);
 //   hessProtonate(atomsLig);
        
//    int idxLig = 0;

 //   res = hessSEAddMol(se_lig, atomsLig, optModeLig, &idxLig);  
    
    // read complex
  //  res = hessSEAddMol(se_complex, atomsProt, optModeProt, &(indices[0]));
   // res = hessSEAddMol(se_complex, atomsLig, optModeLig, &(indices[1]));
    
   // res = hessSERun(se_complex);
    res = hessSERun(se_prot);
   // res = hessSERun(se_lig);
    
    
    double energy_complex, energy_protein, energy_ligand;
    //res = hessSEGetEnergy(se_complex, &energy_complex);
    res = hessSEGetEnergy(se_prot, &energy_protein);
    //res = hessSEGetEnergy(se_lig, &energy_ligand);
   // printf("energy complex: %f\n", energy_complex);
    printf("energy protein: %f\n", energy_protein);


    void *opted_mol;
    res = hessSEGetMol(se_prot, 0, &opted_mol);
    hessSaveMol(opted_mol, params.out[0]);

    
    //printf("energy ligand: %f\n", energy_ligand);
    //printf("binding energy: %f\n", energy_complex - energy_protein - energy_ligand);
    
//    for(int i = 0; i < params.n; i++){
//        void *opted_mol;
//        res = hessSEGetMol(se_complex, indices[i], &opted_mol);
//        hessSaveMol(opted_mol, params.out[i]);
//    }
    
//    printf("Saved molecules with new geoms\n");
    
    free(indices);
    //free(se_complex);
    free(se_prot);
    //free(se_lig);
    free(parser);

    return 0;
}

