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

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "interface/hess.h"

typedef struct {
  char* rec;
  char* lig;
  char* out;
  char* box;
  int depth, number_of_iterations;
  double size_x, size_y, size_z, xc, yc, zc;
  int auto_box;
  int score_only;
  char* optimize;
  int top;
  double granularity;
  unsigned seed;
  int grid_deriv_flag;
} Opt;

int OptCmdParse(int argc, char *argv[], Opt* params) {
  const char *lig_opt = "-l";
  const char *rec_opt = "-r";
  const char *lig_opt_long = "--ligand";
  const char *rec_opt_long = "--receptor";
  const char *out_opt = "-o";
  const char *xc_opt = "--center_x";
  const char *yc_opt = "--center_y";
  const char *zc_opt = "--center_z";
  const char *depth_opt = "--depth";
  const char *number_of_iterations_opt = "--number_of_iterations";
  const char *size_x_opt = "--size_x";
  const char *size_y_opt = "--size_y";
  const char *size_z_opt = "--size_z";
  const char *auto_box_opt = "--autobox_ligand";
  const char *score_only_opt = "--score_only";
  const char *grid_deriv_opt = "--grid_deriv";
  const char *granularity_opt = "--granularity";
  const char *top_opt = "--top";
  const char *optimize_opt = "--optimize";
  const char *seed_opt = "--seed";
  int seed_flag = 0;

  params->rec = NULL;
  params->lig = NULL;
  params->out = NULL;
  params->box = NULL;
  params->optimize = "mc";

  params->grid_deriv_flag = 0;
  params->depth = 10000;
  params->number_of_iterations = 10;
  params->size_x = 40.0;
  params->size_y = 40.0;
  params->size_z = 40.0;
  params->xc = 0.0;
  params->yc = 0.0;
  params->zc = 0.0;
  params->granularity = 0.375;
  params->auto_box = 0;
  params->score_only = 0;
  params->top = 10;

  const char *help_message =
          "Input:\n"
          "-r [ --receptor ] arg         path to the receptor molecule\n"
          "-l [ --ligand ] arg           path to the ligand molecule\n"
          "--autobox_ligand arg          ligand to use for autobox\n"
          "Search space:\n"
          "--center_x arg                X coordinate of the center\n"
          "--center_x arg                Y coordinate of the center\n"
          "--center_x arg                Z coordinate of the center\n"
          "--size_x arg                  size in the X dimension (Angstroms)\n"
          "--size_y arg                  size in the Y dimension (Angstroms)\n"
          "--size_z arg                  size in the Z dimension (Angstroms)\n"
          "Minimization options:\n"
          "--optimize arg                parameter for the global search algorithm ('mc' for Monte Carlo optimisation or 'mc_metropolis' for Monte Carlo with Metropolis criterio)\n"
          "--top arg                     number of top scores\n"
          "--depth arg                   search depth ( number of LBFGS local search runs)\n"
          "--granularity arg             grid splitting coefficient (the smaller the better the approximation of the estimation function)\n"
          "--score_only                  score provided ligand pose\n"
          "--grid_deriv                  use a grid with previous calculation of function gradients\n"
          "Output: \n"
          "-o arg                        output file name, sdf format is used\n"
          "Misc:\n"
          "--seed arg                    explicit random seed\n"
          ;

  for (int i = 0; i < argc; i++) {
    const char* curr = argv[i];
    if (strcmp(curr, "-h") == 0 || strcmp(curr, "--help") == 0) {
      fprintf(stderr, help_message);
      return -1;
    } else if (strcmp(curr, grid_deriv_opt) == 0) {
      params->grid_deriv_flag = 1;
      continue;
    } else if (strcmp(curr, score_only_opt) == 0) {
      params->score_only = 1;
      continue;
    }
    if ((strlen(curr) && curr[0] != '-') || (strlen(curr) > 1 && curr[0] != '-' && curr[1] != '-'))
      continue;
    if ((strcmp(curr, lig_opt) == 0 || strcmp(curr, lig_opt_long) == 0) && i + 1 < argc)
      params->lig = strdup(argv[i + 1]);
    else if ((strcmp(curr, rec_opt) == 0 || strcmp(curr, rec_opt_long) == 0) && i + 1 < argc)
      params->rec = strdup(argv[i + 1]);
    else if (strcmp(curr, out_opt) == 0 && i + 1 < argc)
      params->out = strdup(argv[i + 1]);
    else if (strcmp(curr, xc_opt) == 0 && i + 1 < argc)
      params->xc = atof(argv[i + 1]);
    else if (strcmp(curr, yc_opt) == 0 && i + 1 < argc)
      params->yc = atof(argv[i + 1]);
    else if (strcmp(curr, zc_opt) == 0 && i + 1 < argc)
      params->zc = atof(argv[i + 1]);
    else if (strcmp(curr, depth_opt) == 0 && i + 1 < argc)
      params->depth = atoi(argv[i + 1]);
    else if (strcmp(curr, number_of_iterations_opt) == 0 && i + 1 < argc)
      params->number_of_iterations = atoi(argv[i + 1]);
    else if (strcmp(curr, size_x_opt) == 0 && i + 1 < argc)
      params->size_x = atof(argv[i + 1]);
    else if (strcmp(curr, size_y_opt) == 0 && i + 1 < argc)
      params->size_y = atof(argv[i + 1]);
    else if (strcmp(curr, size_z_opt) == 0 && i + 1 < argc)
      params->size_z = atof(argv[i + 1]);
    else if (strcmp(curr, auto_box_opt) == 0 && i + 1 < argc) {
      params->box = strdup(argv[i + 1]);
      params->auto_box = 1;
    } else if (strcmp(curr, granularity_opt) == 0 && i + 1 < argc)
      params->granularity = atof(argv[i + 1]);
    else if (strcmp(curr, top_opt) == 0 && i + 1 < argc)
      params->top = atoi(argv[i + 1]);
    else if (strcmp(curr, optimize_opt) == 0 && i + 1 < argc)
      params->optimize = strdup(argv[i + 1]);
    else if (strcmp(curr, seed_opt) == 0 && i + 1 < argc) {
      params->seed = atoi(argv[i + 1]);
      seed_flag = 1;
    } else {
      fprintf(stderr,
              "Unrecognized option %s.\n\nCorrect usage:\n\n", help_message);
      return -1;
    }
  }
  if (strcmp(params->optimize, "mc") != 0 && strcmp(params->optimize, "mc_metropolis") != 0) {
    fprintf(stderr, "Unrecognized option %s.'\n\nCorrect usage:\n\n", help_message);
    return -1;
  }
  if (!seed_flag) {
    unsigned time_ui = time(NULL);
    params->seed = time_ui;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  double size_x = 40.0, size_y = 40.0, size_z = 40.0, xc = 0.0, yc = 0.0, zc = 0.0, granularity = 0.375;
  int depth = 10000, number_of_iterations = 100, tops_count = 10, grid_deriv_flag = 0;
  unsigned seed = 12345;
  hessInit();
  hessSetStream(stdout);
  FILE *stream = hessGetStream();
  Opt params;
  if (OptCmdParse(argc, argv, &params) != 0) {
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  size_x = params.size_x;
  size_y = params.size_y;
  size_z = params.size_z;
  xc = params.xc;
  yc = params.yc;
  zc = params.zc;
  depth = params.depth;
  number_of_iterations = params.number_of_iterations;
  const char *rec_path = params.rec;
  const char *lig_path = params.lig;
  const char *out_path = params.out;
  const char *box_path = params.box;
  const char *optimize = params.optimize;
  granularity = params.granularity;
  tops_count = params.top;
  seed = params.seed;
  grid_deriv_flag = params.grid_deriv_flag;
  void* parser = hessCreateParser();
  void* rec_atoms = hessLoadMolecule(parser, rec_path);
  if (rec_atoms == NULL) {
    hessDestroy(parser);
    fprintf(stderr, "Error reading the receptor file. ");
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  void* lig_atoms = hessLoadMolecule(parser, lig_path);
  if (lig_atoms == NULL) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    fprintf(stderr, "Error reading the ligand file. ");
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  void* box_atoms;
  if (box_path != NULL) {
    box_atoms = hessLoadMolecule(parser, box_path);
    if (box_atoms == NULL) {
      hessDestroy(parser);
      hessDestroy(rec_atoms);
      hessDestroy(lig_atoms);
      fprintf(stderr, "Error reading the crystal ligand file. ");
      fprintf(stderr, "%s\n", hessGetLastError());
      return -1;
    }
  }
  double box[6] = {xc, yc, zc, size_x, size_y, size_z};
  if (params.auto_box) {
    hessCalcAutobox(box_atoms, box);
  }
  if (hessDeleteHydrogens(rec_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessDeleteHydrogens(lig_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessDeleteFormalCharges(rec_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessDeleteFormalCharges(lig_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  hessSetPh(rec_atoms, 7.4);
  hessSetPh(lig_atoms, 7.4);
  if (hessAssignBonds(rec_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessAssignBonds(lig_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessProtonate(rec_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessProtonate(lig_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessAddHydrogens(rec_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (hessAddHydrogens(lig_atoms) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  double result[2] = {0};
  fprintf(stream, "Version: Dec 19, 2023\n");
  void* opt_molecule = hessMakeOptimizableMolecule(lig_atoms, rec_atoms, box, optimize, granularity, seed);
  if (opt_molecule == NULL) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  if (params.score_only) {
    hessWriteScoreOnly(opt_molecule);
    hessDestroy(lig_atoms);
    hessDestroy(rec_atoms);
    hessDestroy(box_atoms);
    return 0;
  }
  clock_t grid_time, opt_time;
  grid_time = clock();
  hessFillGrid(opt_molecule, grid_deriv_flag);
  grid_time = clock() - grid_time;
  fprintf(stream, "Grid built in %.2f s\n", ((double) grid_time) / CLOCKS_PER_SEC);

  opt_time = clock();
  if (hessRunOptimize(number_of_iterations, depth, opt_molecule, result, tops_count) != 0) {
    hessDestroy(parser);
    hessDestroy(rec_atoms);
    hessDestroy(lig_atoms);
    hessDestroy(box_atoms);
    hessDestroy(opt_molecule);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  opt_time = clock() - opt_time;
  fprintf(stream, "Optimization took %.2f seconds\n", ((double) opt_time) / CLOCKS_PER_SEC);
  fprintf(stream, "Best intermolecular and total binding energies:\n");
  fprintf(stream, "%7.3f %7.3f\n", result[0], result[0] + result[1]);
  if (out_path != NULL && hessSaveSdf(opt_molecule, out_path) != 0) {
    hessDestroy(parser);
    hessDestroy(lig_atoms);
    hessDestroy(rec_atoms);
    hessDestroy(box_atoms);
    hessDestroy(opt_molecule);
    fprintf(stderr, "%s\n", hessGetLastError());
    return -1;
  }
  hessDestroy(lig_atoms);
  hessDestroy(rec_atoms);
  hessDestroy(opt_molecule);
  hessDestroy(box_atoms);
  hessDestroy(parser);
  return 0;
}
