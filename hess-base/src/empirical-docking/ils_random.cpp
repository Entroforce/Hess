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

#include "global_optimizer.h"
#include "optimization_methods.h"
#include "LBFGS.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXf;
using namespace LBFGSpp;

Eigen::VectorXd ils_random(Optimizable_molecule& mol, int depth, double dif) {
  int solve_size = mol.encoding.size();
  double size_x = mol.size_x;
  double size_y = mol.size_y;
  double size_z = mol.size_z;
  double conr[3] = {size_x, size_y, size_z};
  double start_conr[3] = {size_x, size_y, size_z};
  int conr_id = 0;
  double fx = 1000;
  VectorXd x = VectorXd::Zero(solve_size);
  VectorXd bestans = VectorXd::Zero(solve_size);
  int start_pos_attempts = 0;
  int MAX_START_ATTEMPTS = 40;
  int MAX_START_SECOND_ATTEMPTS = 80;
  bool (*box_start_checker)(hess::Molecule *lig, double size_x, double size_y, double size_z, simplified_tree& tr, int& ex_count, const vector<int>& encoding_inv, const Eigen::VectorXd & x);
  box_start_checker = check_exceeded_box_limits_start;
  bool was_moved = false;
  int ex_count = 0;
  while (true) {
    conr_id = 0;
    start_pos_attempts += 1;
    if (start_pos_attempts > MAX_START_ATTEMPTS && !was_moved) {
      start_conr[0] = size_x / 2;
      start_conr[1] = size_y / 2;
      start_conr[2] = size_z / 2;
      box_start_checker = check_exceeded_box_limits;
      was_moved = true;
    }
    if (start_pos_attempts > MAX_START_SECOND_ATTEMPTS && was_moved) {
      throw HessException("Failed to place the ligand in the box. Check box dimensions. Maybe he's too small");
    }
    for (int i = 0; i < solve_size - 3; i++)
      x[i] = random_angle_ils();
    for (int i = solve_size - 3; i < solve_size; i++) {
      x[i] = random_number_ils(start_conr[conr_id]);
      conr_id++;
    }
    if (!box_start_checker(mol.ligand, size_x, size_y, size_z, mol.tr, ex_count, mol.encoding_inv, x))
      break;
  }
  LBFGSParam<double> param;
  LBFGSSolver<double> solver(param);
  param.epsilon = 0.0001;
  param.epsilon_rel = 0.0001;
  double alpha = 0.99995, beta = 0.99995;
  double best = 1000.0;
  for (int i = 0; i < depth; i++) {
    if (fx < best) {
      bestans.noalias() = x;
      best = min(best, fx);
    }
    VectorXd new_x1 = VectorXd::Zero(x.size());
    random_change(new_x1, x, dif, conr);
    solver.minimize(mol, new_x1, fx);
    x.noalias() = new_x1;
    param.epsilon *= alpha;
    param.epsilon_rel *= alpha;
    dif *= beta;
  }
  return bestans;
}