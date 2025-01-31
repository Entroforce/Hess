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

#ifndef EIGEN
#define EIGEN
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#endif


#include "scoring_function/features.h"
#include "model/simplified_tree.h"
#include "model/array3d.h"

using Eigen::VectorXd;

class Optimizable_molecule {
public:
  Optimizable_molecule(hess::Molecule* _lig, hess::Molecule* _prot, double* box, const char* _optimize, double _granularity, unsigned _seed);
  void fill_grid();
  void fill_grid_with_derivs();
  double calc_derivative_fast_type_0_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
  double calc_derivative_fast_type_1_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
  double calc_derivative_fast_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
  void set_max_bfgs_iterations(int max_bfgs_iterations);
  void assign_hunt_cap();
  void assign_aut_cap();
  double operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad);
  ~Optimizable_molecule();


  hess::Molecule *ligand;
  hess::Molecule *receptor;
  simplified_tree tr;
  VectorXd fixedcon;
  int type;
  vector<int> encoding;
  vector<int> encoding_inv;
  double epsl = 0.0001;
  int grid_step = 0;
  double size_x;
  double size_y;
  double size_z;
  double xc;
  double yc;
  double zc;
  double granularity = 0.375;
  vector<int> rotsmap;
  vector<Array3D<double>> grid_new;
  vector<Array3D<hess::Vec3d>> grid_deriv_new;
  double gd[3][3] = {
    {0},};
  double init[3] = {0};
  double corner1[3] = {0};
  double corner2[3] = {0};
  hess::Vec3d init_vec;
  double range[3] = {0};
  double factor[3] = {0};
  hess::Vec3d factor_vec;
  double dim_fl_minus_1[3] = {0};
  double factor_inv[3] = {0};
  hess::Vec3d factor_inv_vec;
  int dim_x = 0, dim_y = 0, dim_z = 0;
  set<smt> ligand_types;
  std::vector<smt> ligand_types_vec;
  double slope = 1e6;
  ConfIndependentInputs in;
  int rot_bonds_count = 0;
  int num_steps = 0;
  int mut_amplitude = 2;
  int temperature = 1.2;
  int num_saved_mins = 20;
  double min_rmsd = 1.0;
  int hunt_cap[2] = {10, 10};
  int aut_cap[2] = {1000, 1000};
  int cap[2] = {1000, 1000};
  int solve_count = 0;
  unsigned max_bfgs_iterations;
  const char* scoring = "vinardo";
  bool grid_with_derivs = false;
  int tops_count = 10;
  string optimize = "mc";
  vector<hess::Molecule*> result_mols;
  unsigned seed;
  double ligand_center[3] = {0};
};

void find_ligand_pairs(hess::Molecule* lig);

