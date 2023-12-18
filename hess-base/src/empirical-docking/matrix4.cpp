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

#include <string.h>

#include "matrix4.h"
using std::array;

Matrix4::Matrix4() {
}

Matrix4::~Matrix4() {
}

Matrix4::Matrix4(const Matrix4& o) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      this->matr[i][j] = o.matr[i][j];
}

Matrix4::Matrix4(double elem[]) {
  memcpy(this->matr, elem, sizeof (this->matr));
}

Matrix4& Matrix4::operator=(const Matrix4& o) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      this->matr[i][j] = o.matr[i][j];
  return *this;
}

Matrix4& Matrix4::operator=(Matrix4&& o) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      this->matr[i][j] = o.matr[i][j];
      o.matr[i][j] = 0;
    }
  return *this;
}

Matrix4::Matrix4(Matrix4&& o) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      this->matr[i][j] = o.matr[i][j];
      o.matr[i][j] = 0;
    }
}

Matrix4 mlt(const Matrix4& a, const Matrix4& b) {
  Matrix4 nw;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      for (int h = 0; h < 4; h++)
        nw.matr[i][h] += a.matr[i][j] * b.matr[j][h];
  return nw;
}

Matrix4 transpose(Matrix4& a) {
  Matrix4 nw;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      nw.matr[i][j] = a.matr[j][i];
  return nw;
}

Matrix4 shift_by_vec(double x, double y, double z) {
  Matrix4 res;
  res.matr[0][0] = 1.0;
  res.matr[3][0] = x;
  res.matr[1][1] = 1.0;
  res.matr[3][1] = y;
  res.matr[2][2] = 1.0;
  res.matr[3][2] = z;
  res.matr[3][3] = 1.0;

  return res;
}

array<double, 3>vecprod(array<double, 3>& v1, array<double, 3>& v2) {
  array<double, 3>res;
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  return res;
}

void normalize(array<double, 3>& v) {
  double norm = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  norm = sqrt(norm);
  v[0] /= norm;
  v[1] /= norm;
  v[2] /= norm;
}

Matrix4 getbasis(double x, double y, double z) {
  array<double, 3> p1 = {x, y, z}, p2 = {0.1, 0.2, 0.3};
  auto p3 = vecprod(p1, p2);
  p2 = vecprod(p1, p3);
  normalize(p1);
  normalize(p2);
  normalize(p3);
  Matrix4 res;
  res.matr[3][3] = 1;
  res.matr[0][0] = p1[0];
  res.matr[0][1] = p1[1];
  res.matr[0][2] = p1[2];
  res.matr[1][0] = p2[0];
  res.matr[1][1] = p2[1];
  res.matr[1][2] = p2[2];
  res.matr[2][0] = p3[0];
  res.matr[2][1] = p3[1];
  res.matr[2][2] = p3[2];
  return res;
}

Matrix4 rotate_by_x_vec(double angle) {
  Matrix4 res;
  res.matr[0][0] = 1.0;
  res.matr[3][3] = 1.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[1][1] = cs, res.matr[1][2] = -sn, res.matr[2][1] = sn, res.matr[2][2] = cs;
  return res;
}

Matrix4 global_rotate_by_x(double angle) {
  Matrix4 res;
  res.matr[0][0] = 1.0;
  res.matr[3][3] = 1.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[1][1] = cs, res.matr[1][2] = -sn, res.matr[2][1] = sn, res.matr[2][2] = cs;
  return res;
}

Matrix4 global_rotate_by_y(double angle) {
  Matrix4 res;
  res.matr[1][1] = 1.0;
  res.matr[3][3] = 1.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[0][0] = cs, res.matr[0][2] = sn, res.matr[2][0] = -sn, res.matr[2][2] = cs;
  return res;
}

Matrix4 global_rotate_by_z(double angle) {
  Matrix4 res;
  res.matr[3][3] = 1.0;
  res.matr[2][2] = 1.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[0][0] = cs, res.matr[0][1] = -sn, res.matr[1][0] = sn, res.matr[1][1] = cs;
  return res;
}

Matrix4 global_rotate_by_xyz(double a, double b, double g) {
  Matrix4 res;
  double sn_a = sin(a), cs_a = cos(a);
  double sn_b = sin(b), cs_b = cos(b);
  double sn_g = sin(g), cs_g = cos(g);
  res.matr[0][0] = cs_b * cs_g, res.matr[0][1] = -sn_g * cs_b, res.matr[0][2] = sn_b;
  res.matr[1][0] = sn_a * sn_b * cs_g + sn_g * cs_a, res.matr[1][1] = -sn_a * sn_b * sn_g + cs_a * cs_g, res.matr[1][2] = -sn_a * cs_b;
  res.matr[2][0] = sn_a * sn_g - sn_b * cs_a * cs_g, res.matr[2][1] = sn_a * cs_g + sn_b * sn_g * cs_a, res.matr[2][2] = cs_a * cs_b;
  return res;
}

Matrix4 global_rotate_by_x_deriv(double angle) {
  Matrix4 res;
  res.matr[0][0] = 0.0;
  res.matr[3][3] = 0.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[1][1] = -sn, res.matr[1][2] = -cs, res.matr[2][1] = cs, res.matr[2][2] = -sn;
  return res;
}

Matrix4 global_rotate_by_y_deriv(double angle) {
  Matrix4 res;
  res.matr[1][1] = 0.0;
  res.matr[3][3] = 0.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[0][0] = -sn, res.matr[0][2] = cs, res.matr[2][0] = -cs, res.matr[2][2] = -sn;
  return res;
}

Matrix4 global_rotate_by_z_deriv(double angle) {
  Matrix4 res;
  res.matr[3][3] = 0.0;
  res.matr[2][2] = 0.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[0][0] = -sn, res.matr[0][1] = -cs, res.matr[1][0] = cs, res.matr[1][1] = -sn;
  return res;
}

Matrix4 rotate(double x, double y, double z, double angle) {
  auto basis = getbasis(x, y, z);
  Matrix4 res = transpose(basis);
  res = mlt(res, rotate_by_x_vec(angle));
  res = mlt(res, basis);
  return res;
}

void rotate_matrix_to_euler_ZYX(Matrix4& M_new, double* angles) {
  double r00 = M_new.matr[0][0];
  double r01 = M_new.matr[0][1];
  double r02 = M_new.matr[0][2];
  double r10 = M_new.matr[1][0];
  double r11 = M_new.matr[1][1];
  double r12 = M_new.matr[1][2];
  double r20 = M_new.matr[2][0];
  double r21 = M_new.matr[2][1];
  double r22 = M_new.matr[2][2];
  double x_angle = 0;
  double y_angle = 0;
  double z_angle = 0;
  if (r20 < 1) {
    if (r20 > -1) {
      y_angle = asin(-r20);
      z_angle = atan2(r10, r00);
      x_angle = atan2(r21, r22);
    } else {
      y_angle = M_PI / 2;
      z_angle = -atan2(-r12, r11);
      x_angle = 0;
    }
  } else {
    y_angle = -M_PI / 2;
    z_angle = atan2(-r12, r11);
    x_angle = 0;
  }
  angles[0] = x_angle;
  angles[1] = y_angle;
  angles[2] = z_angle;
}
