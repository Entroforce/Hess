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

#pragma once

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  void hessInit();
  int hessAssignBonds(void *old_mol);
  int hessProtonate(void *mol);
  int hessAddHydrogens(void* mol);
  void* hessCreateParser();
  void* hessLoadMolecule(void* parser, const char* file_name);
  int hessSaveMol(void* mol, const char* path);
  int hessSaveSdf(void* mols_void, const char *path);
  void hessDestroy(void* object);
  const char* hessGetLastError();
  void hessSetPh(void *mol, double pH);
  void hessSetLogFile(const char * path);
  FILE* hessGetStream();
  void hessSetStream(FILE* stream);
  int hessDeleteHydrogens(void *mol);
  int hessDeleteFormalCharges(void *mol);

  void hessCalcAutobox(void* ligand_v, double* box);
  void hessWriteScoreOnly(void* opt_mol, double* receptor_center);
  int hessRunRandomIls(int number_of_iterations, int depth, void* opt_v, double* result_array, int tops_count);
  int hessRunOptimize(int number_of_iterations, int depth, void* opt_v, double* result_array, int tops_count);
  int hessRunMonteCarlo(int number_of_iterations, void* opt_v, double* result_array, int tops_count);
  void* hessMakeOptimizableMolecule(void* ligand_v, void* rec_v, double* box, const char* optimize, double granularity, unsigned seed);
  void hessFillGridWithDerivs(void* opt_v);
  void hessFillGrid(void* opt_v, int grid_deriv_flag);

#ifdef __cplusplus
}
#endif

