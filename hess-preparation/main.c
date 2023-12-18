/**************************************************************************
 * This file is part of the Hess project
 * Copyright (C) 2023 Entroforce LLC
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

#include "interface/hess.h"

int main(int argc, char *argv[]) {

  hessInit();

  if (argc != 3) {
    fprintf(stderr, "Usage: hess-preparation <input pdb file> <output molfile>\n");
    return -1;
  }
  const char *path = argv[1];
  void *parser = hessCreateParser();

  void *atoms = hessLoadMolecule(parser, path);
  if (atoms == NULL) {
    hessDestroy(parser);
    fprintf(stderr, "Error in reading the molecule file. ");
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }

  hessSetPh(atoms, 7.4);

  if (hessDeleteHydrogens(atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }

  if (hessDeleteFormalCharges(atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }

  if (hessAssignBonds(atoms) != 0) {
    fprintf(stderr, "%s\n", hessGetLastError());
    hessDestroy(parser);
    hessDestroy(atoms);
    return -1;
  }
  if (hessProtonate(atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessAddHydrogens(atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessSaveMol(atoms, argv[2]) != 0) {
    fprintf(stderr, "%s\n", hessGetLastError());
    hessDestroy(parser);
    hessDestroy(atoms);
    return -1;
  }

  hessDestroy(parser);
  hessDestroy(atoms);
  return 0;
}
