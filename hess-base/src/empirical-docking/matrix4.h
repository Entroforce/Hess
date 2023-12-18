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

#include <math.h>
#include "molecule.h"

class Matrix4 {
public:
  double matr[4][4] = {0};
  Matrix4();
  Matrix4(double elem[]);
  Matrix4(const Matrix4& o);
  Matrix4& operator=(const Matrix4& o);
  Matrix4& operator=(Matrix4&& o);
  Matrix4(Matrix4&& o);
  ~Matrix4();

  inline static const Matrix4 & identity() {
    static double _arr[] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    static Matrix4 _instance(_arr);
    return _instance;
  }

};


Matrix4 transpose(Matrix4& a);
Matrix4 getbasis(double x, double y, double z);
Matrix4 mlt(const Matrix4& a, const Matrix4& b);
Matrix4 shift_by_vec(double x, double y, double z);
Matrix4 rotate(double x, double y, double z, double angle);
Matrix4 global_rotate_by_x(double angle);
Matrix4 global_rotate_by_y(double angle);
Matrix4 global_rotate_by_z(double angle);
Matrix4 global_rotate_by_x_deriv(double angle);
Matrix4 global_rotate_by_y_deriv(double angle);
Matrix4 global_rotate_by_z_deriv(double angle);
Matrix4 global_rotate_by_xyz(double a, double b, double g);
void rotate_matrix_to_euler_ZYX(Matrix4& M_new, double* angles);
