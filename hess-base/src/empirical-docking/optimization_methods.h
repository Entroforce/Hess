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

#include "transformation.h"
#include "array3d.h"
using Eigen::VectorXd;

void calc_all_da_dalpha(vector<hess::Vec3d> &da_dalpha, hess::Molecule *lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x, vector<int>& rotsmap);
void calc_all_da_dalpha_type1(vector<hess::Vec3d>& da_dalpha, hess::Molecule *lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x);

int grid_encode(int a, int x, int y, int z, const int sz);
int encode(int a, int x, int y, int z, const int sz);
double random_number_ils(double q);
void random_change_type_0(VectorXd& v, VectorXd& old_v, double coef, double conr[3]);
void random_change_type_2(VectorXd& v, VectorXd& old_v, double coef, double conr[3]);
void random_change(VectorXd& v, VectorXd& old_v, double coef, double conr[3]);
double random_angle_ils();
