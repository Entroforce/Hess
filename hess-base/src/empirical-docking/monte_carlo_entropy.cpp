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
#include "LBFGS.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXf;
using namespace LBFGSpp;

void monte_carlo(Optimizable_molecule& mol, vector<mc_out>& outs) {
  LBFGSParam<double> param;
  LBFGSSolver<double> solver(param);
  param.epsilon = 0.0001;
  param.epsilon_rel = 0.0001;
  param.max_iterations = mol.max_bfgs_iterations;
  int solve_size = mol.rot_bonds_count + 6;
  VectorXd tmp = VectorXd::Zero(solve_size);
  VectorXd candidate = VectorXd::Zero(solve_size);
  first_randomize_monte_carlo(tmp, mol);
  double best_e = max_fl;
  double tmp_e = 1000;
  double can_e = 1000;
  int num_steps = mol.num_steps;
  for (int step = 0; step < num_steps; step++) {
    candidate.noalias() = tmp;
    mutate_solve(candidate, mol);
    mol.assign_hunt_cap();
    solver.minimize(mol, candidate, can_e);
    if (step == 0 || metropolis_accept(tmp_e, can_e, mol.temperature)) {
      tmp.noalias() = candidate;
      tmp_e = can_e;
      if (tmp_e < best_e || outs.size() < mol.num_saved_mins) {
        mol.assign_aut_cap();
        solver.minimize(mol, tmp, tmp_e);
        add_to_output_container(outs, tmp, mol, tmp_e); // 20 - max size
        if (tmp_e < best_e)
          best_e = tmp_e;
      }
    }
  }
}
