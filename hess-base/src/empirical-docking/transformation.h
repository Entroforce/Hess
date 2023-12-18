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
#include "molecule.h"
#include "matrix4.h"
#include "simplified_tree.h"

#ifndef EIGEN
#define EIGEN
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#endif

void tree_init(simplified_tree &tr, hess::Molecule *pd);
void transformation(hess::Molecule *pd, simplified_tree &tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x);
void rec(hess::Molecule *pd, simplified_tree &tr, int v, Matrix4 curmatr, const vector<int>& encoding_inv, const Eigen::VectorXd& x);
Matrix4 rotate_by_x_vec_deriv(double angle);
Matrix4 rotate_deriv(double x, double y, double z, double angle);
void rec_deriv(hess::Molecule *pd, simplified_tree &tr, vector<double>&rotations, int v, Matrix4 curmatr, Matrix4 curder, int alpha_id, int meet_alpha_id, vector<hess::Vec3d>& deriva);
