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

#include "fast_gradient.h"
#include "ils_primitives.h"
#include "simplified_tree.h"

vector <double> calc_autobox(hess::Molecule *ligand);
void set_insubtree(hess::Molecule* lig);
void calc_center_mass_frag(hess::Molecule *pd, simplified_tree &tr, int v, int& count, hess::Vec3d& average);
void set_parents_to_tree(hess::Molecule* lig);
void calc_global_shifts(vector<vector<double>>&global_shifts, Optimizable_molecule& mol, int niter);
void form_ils_results(vector<pair<Eigen::VectorXd, pair<double, double>>>& result_pairs_sort, Optimizable_molecule* opt);
pair<double, double> calc_energy_for_result(hess::Molecule *lig, hess::Molecule *prot, simplified_tree &tr, ConfIndependentInputs& in, const vector<int>& encoding_inv, const Eigen::VectorXd& x, const char* scoring);
double calc_energy(hess::Molecule *lig, hess::Molecule *prod, ConfIndependentInputs& in, const char* scoring);
void sort_configurations(vector<pair<Eigen::VectorXd, pair<double, double>>>& result_pairs);
