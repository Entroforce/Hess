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

#include "optimizable_molecule.h"

bool check_conf(hess::Molecule *lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x);
bool check_exceeded_box_limits(hess::Molecule *lig, double size_x, double size_y, double size_z, simplified_tree& tr, int& ex_count, const vector<int>& encoding_inv, const Eigen::VectorXd& x);
bool check_exceeded_box_limits_start(hess::Molecule *lig, double size_x, double size_y, double size_z, simplified_tree& tr, int& ex_count, const vector<int>& encoding_inv, const Eigen::VectorXd& x);
Eigen::VectorXd ils_random(Optimizable_molecule& mol, int depth, double dif);

class mc_out {
public:

  mc_out(const vector<hess::Vec3d>& _coords, const Eigen::VectorXd& _solve, double _e = 0.0) : e(_e) {
    solve = _solve;
    coords = _coords;
  }

  mc_out(double _e = 0.0) : e(_e) {
    solve = VectorXd::Zero(0);
    coords = {};
  }
  vector<hess::Vec3d> coords;
  double e;
  Eigen::VectorXd solve;

  double get_energy()const {
    return this->e;
  }

  const vector<hess::Vec3d>& get_coords()const {
    return this->coords;
  }

  const Eigen::VectorXd& get_solve()const {
    return this->solve;
  }
};

void monte_carlo(Optimizable_molecule& mol, vector<mc_out>& outs);
void first_randomize_monte_carlo(Eigen::VectorXd& x, Optimizable_molecule& mol);
void mutate_solve(Eigen::VectorXd& x, Optimizable_molecule& mol);
bool metropolis_accept(double old_f, double new_f, double temperature);
void add_to_output_container(vector<mc_out>& outs, VectorXd& tmp, Optimizable_molecule& mol, double e);
void merge_output_containers(vector<mc_out>& all_outs, vector<mc_out>& outs, Optimizable_molecule& mol);
int random_int_ab(int min, int max);
