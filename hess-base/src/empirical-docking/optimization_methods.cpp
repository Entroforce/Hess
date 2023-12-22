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

#include "optimization_methods.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXf;

const double eps = 0.00001;

double circlecome_ils(double x, double y, double z) {
  double val = 0;
  if (y > x - eps) {
    if (y - x < M_PI)
      val = x * (1 - z) + y * z;
    else {
      y -= 2 * M_PI;
      val = x * (1 - z) + y*z;
      if (val < 0)
        val += 2 * M_PI;
    }
    return val;
  } else {
    return circlecome_ils(y, x, 1 - z);
  }
}

double usualcome_ils(double x, double y, double z) {
  return x * (1 - z) + y*z;
}

double random_angle_ils() {
  double val = (rand() % 10000);
  val /= 10000.0;
  return val * 2.0 * M_PI;
}

double random_number_ils(double q) {
  double val = (rand() % 10000);
  val /= 10000.0;
  return val * q - q / 2;
}

void random_change(VectorXd& v, VectorXd& old_v, double coef, double conr[3]) {
  int conr_id = 0;
  for (int i = 0; i < old_v.size() - 3; i++) {
    double a = random_angle_ils();
    v[i] = circlecome_ils(old_v[i], a, coef);
  }
  for (int i = v.size() - 3; i < v.size(); i++) {
    double a = random_number_ils(conr[conr_id]);
    v[i] = usualcome_ils(old_v[i], a, coef);
    conr_id++;
  }
}

void calc_all_da_dalpha_type1(vector<hess::Vec3d> &da_dalpha, hess::Molecule *lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  int trs_count = lig->get_rotatable_bonds_count();
  rec(lig, tr, tr.rootid, Matrix4::identity(), encoding_inv, x);
  Matrix4 M1 = global_rotate_by_x(x[trs_count]);
  Matrix4 DM1 = global_rotate_by_x_deriv(x[trs_count]);
  Matrix4 M2 = global_rotate_by_y(x[trs_count + 1]);
  Matrix4 DM2 = global_rotate_by_y_deriv(x[trs_count + 1]);
  Matrix4 M3 = global_rotate_by_z(x[trs_count + 2]);
  Matrix4 DM3 = global_rotate_by_z_deriv(x[trs_count + 2]);
  Matrix4 D[3];
  Matrix4 MALL = mlt(M1, M2);
  MALL = mlt(MALL, M3);
  D[0] = mlt(DM1, M2);
  D[0] = mlt(D[0], M3);
  D[1] = mlt(M1, DM2);
  D[1] = mlt(D[1], M3);
  D[2] = mlt(M1, M2);
  D[2] = mlt(D[2], DM3);
  for (int i = 0; i < 3; i++) {
    const Matrix4 &trf = D[i];
    for (int j = 0; j < tr.sz; j++) {
      da_dalpha[tr.sz * i + j].x = lig->atoms[j].x * trf.matr[0][0] + lig->atoms[j].y * trf.matr[1][0] + lig->atoms[j].z * trf.matr[2][0] + trf.matr[3][0];
      da_dalpha[tr.sz * i + j].y = lig->atoms[j].x * trf.matr[0][1] + lig->atoms[j].y * trf.matr[1][1] + lig->atoms[j].z * trf.matr[2][1] + trf.matr[3][1];
      da_dalpha[tr.sz * i + j].z = lig->atoms[j].x * trf.matr[0][2] + lig->atoms[j].y * trf.matr[1][2] + lig->atoms[j].z * trf.matr[2][2] + trf.matr[3][2];
    }
  }
  for (int i = 0; i < tr.sz; i++) {
    double newx = lig->atoms[i].x * MALL.matr[0][0] + lig->atoms[i].y * MALL.matr[1][0] + lig->atoms[i].z * MALL.matr[2][0] + x[encoding_inv[tr.sz + 3]];
    double newy = lig->atoms[i].x * MALL.matr[0][1] + lig->atoms[i].y * MALL.matr[1][1] + lig->atoms[i].z * MALL.matr[2][1] + x[encoding_inv[tr.sz + 4]];
    double newz = lig->atoms[i].x * MALL.matr[0][2] + lig->atoms[i].y * MALL.matr[1][2] + lig->atoms[i].z * MALL.matr[2][2] + x[encoding_inv[tr.sz + 5]];
    lig->atoms[i].x = newx;
    lig->atoms[i].y = newy;
    lig->atoms[i].z = newz;
  }
}

void calc_all_da_dalpha(vector<hess::Vec3d>& da_dalpha, hess::Molecule* lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x, vector<int>& rotsmap) {
  int rotsmap_size = rotsmap.size();
  rec(lig, tr, tr.rootid, Matrix4::identity(), encoding_inv, x);
  int trs_count = rotsmap_size;
  Matrix4 M1 = global_rotate_by_x(x[trs_count]);
  Matrix4 DM1 = global_rotate_by_x_deriv(x[trs_count]);
  Matrix4 M2 = global_rotate_by_y(x[trs_count + 1]);
  Matrix4 DM2 = global_rotate_by_y_deriv(x[trs_count + 1]);
  Matrix4 M3 = global_rotate_by_z(x[trs_count + 2]);
  Matrix4 DM3 = global_rotate_by_z_deriv(x[trs_count + 2]);
  Matrix4 D[3];
  Matrix4 MALL(mlt(M1, M2));
  MALL = mlt(MALL, M3);
  D[0] = mlt(DM1, M2);
  D[0] = mlt(D[0], M3);
  D[1] = mlt(M1, DM2);
  D[1] = mlt(D[1], M3);
  D[2] = mlt(M1, M2);
  D[2] = mlt(D[2], DM3);
  const vector<vector<int>> &insubt = lig->get_insubt();
  const vector<int> &tree_parents = lig->get_tree_parents();
  for (int i = 0; i < rotsmap_size; i++) {
    int pos = rotsmap[i];
    int father = tree_parents[i];
    Matrix4 trf = shift_by_vec(-lig->atoms[father].x, -lig->atoms[father].y, -lig->atoms[father].z);
    trf = mlt(trf, rotate_deriv(lig->atoms[pos].x - lig->atoms[father].x, lig->atoms[pos].y - lig->atoms[father].y, lig->atoms[pos].z - lig->atoms[father].z, 0));
    trf = mlt(trf, shift_by_vec(lig->atoms[father].x, lig->atoms[father].y, lig->atoms[father].z));
    trf = mlt(trf, MALL);
    for (int j = 0; j < tr.sz; j++) {
      if (!insubt[pos][j]) {
        da_dalpha[tr.sz * i + j].zero();
      } else {
        da_dalpha[tr.sz * i + j].x = lig->atoms[j].x * trf.matr[0][0] + lig->atoms[j].y * trf.matr[1][0] + lig->atoms[j].z * trf.matr[2][0] + trf.matr[3][0];
        da_dalpha[tr.sz * i + j].y = lig->atoms[j].x * trf.matr[0][1] + lig->atoms[j].y * trf.matr[1][1] + lig->atoms[j].z * trf.matr[2][1] + trf.matr[3][1];
        da_dalpha[tr.sz * i + j].z = lig->atoms[j].x * trf.matr[0][2] + lig->atoms[j].y * trf.matr[1][2] + lig->atoms[j].z * trf.matr[2][2] + trf.matr[3][2];
      }
    }
  }
  for (int i = 0; i < 3; i++) {
    const Matrix4 &trf = D[i];
    for (int j = 0; j < tr.sz; j++) {
      da_dalpha[tr.sz * (rotsmap_size + i) + j].x = lig->atoms[j].x * trf.matr[0][0] + lig->atoms[j].y * trf.matr[1][0] + lig->atoms[j].z * trf.matr[2][0] + trf.matr[3][0];
      da_dalpha[tr.sz * (rotsmap_size + i) + j].y = lig->atoms[j].x * trf.matr[0][1] + lig->atoms[j].y * trf.matr[1][1] + lig->atoms[j].z * trf.matr[2][1] + trf.matr[3][1];
      da_dalpha[tr.sz * (rotsmap_size + i) + j].z = lig->atoms[j].x * trf.matr[0][2] + lig->atoms[j].y * trf.matr[1][2] + lig->atoms[j].z * trf.matr[2][2] + trf.matr[3][2];
    }
  }
  for (int i = 0; i < tr.sz; i++) {
    double newx = lig->atoms[i].x * MALL.matr[0][0] + lig->atoms[i].y * MALL.matr[1][0] + lig->atoms[i].z * MALL.matr[2][0] + x[encoding_inv[tr.sz + 3]];
    double newy = lig->atoms[i].x * MALL.matr[0][1] + lig->atoms[i].y * MALL.matr[1][1] + lig->atoms[i].z * MALL.matr[2][1] + x[encoding_inv[tr.sz + 4]];
    double newz = lig->atoms[i].x * MALL.matr[0][2] + lig->atoms[i].y * MALL.matr[1][2] + lig->atoms[i].z * MALL.matr[2][2] + x[encoding_inv[tr.sz + 5]];
    lig->atoms[i].x = newx;
    lig->atoms[i].y = newy;
    lig->atoms[i].z = newz;
  }
}