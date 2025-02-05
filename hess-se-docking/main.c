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
    const char *prot_opt = "-popt";
    const char *prot_out = "-pout";

  params->n = 1;
  params->mol = (char **) malloc(params->n * sizeof(char *));
  params->out = (char **) malloc(params->n * sizeof(char *));
  params->optimize = (char **) malloc(params->n * sizeof(char *));

  const char *help_message =
          "Input:\n"
          "-p arg         path to the molecule, pdb format is used\n"
          "Minimization options:\n"
          "-popt arg      parameter for the molecule optimization space ('n' for no optimization; 'y' for all molecule optimization; 'h' for only hydrogens optimization)\n"
          "Output: \n"
          "-pout arg      output molecule file name, sdf format is used\n"
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
      
    else if (strcmp(curr, prot_out) == 0 && i + 1 < argc){
        params->out[0] = strdup(argv[i + 1]);
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
    
    void * se_mol = hessSECreate();
    
    int res;
    
    int optModeMol = 0;
    
    // read mol
    
    if(strcmp(params.optimize[0], "y") == 0) {
        optModeMol = 1;
    } else if(strcmp(params.optimize[0], "h") == 0) {
        optModeMol = 2;
    }
    
    void* atomsMol = hessLoadMolecule(parser, params.mol[0]);
    if (atomsMol == NULL) {
        hessDestroy(parser);
        fprintf(stderr, "Error reading the receptor file. ");
        fprintf(stderr, "%s\n", hessGetLastError());
        return -1;
    }

    hessAssignBonds(atomsMol);
    hessProtonate(atomsMol);
    
    int idxMol = 0;    
    res = hessSEAddMol(se_mol, atomsMol, optModeMol, &idxMol);
    res = hessSERun(se_mol);

    double energy_mol;
    res = hessSEGetEnergy(se_mol, &energy_mol);
    printf("mol energy: %f\n", energy_mol);

    void *opted_mol;
    res = hessSEGetMol(se_mol, 0, &opted_mol);
    hessSaveMol(opted_mol, params.out[0]);

    hessDestroy(parser);
    hessDestroy(atomsMol);

    return 0;
}

