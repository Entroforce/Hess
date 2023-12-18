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

#include "transformation.h"
using namespace std;

simplified_tree::simplified_tree(int n) {
  rootid = 0;
  sz = n;
  children.resize(n);
  bond_to_parent.resize(n);
}

simplified_tree::~simplified_tree() {
  children.clear();
  bond_to_parent.clear();
}

void dfs(simplified_tree &tr, hess::Molecule *pd, vector<int>& used, int v) {
  used[v] = 1;
  for (int i = pd->edgeBegin(); i != pd->edgeEnd(); i = pd->edgeNext(i)) {
    int fir = pd->bond_atom_beg(i);
    int sec = pd->bond_atom_end(i);
    int type = pd->get_bond(i)->rotatable;
    if (fir == v) {
      if (!used[sec]) {
        tr.bond_to_parent[sec] = type;
        tr.children[v].push_back(sec);
        dfs(tr, pd, used, sec);
      }
    }
    if (sec == v) {
      if (!used[fir]) {
        tr.bond_to_parent[fir] = type;
        tr.children[v].push_back(fir);
        dfs(tr, pd, used, fir);
      }
    }
  }
}

void tree_init(simplified_tree &tr, hess::Molecule *pd) {
  vector<int>used(tr.sz);
  used[0] = 1;
  dfs(tr, pd, used, 0);
}

void rec(hess::Molecule *pd, simplified_tree &tr, int v, Matrix4 curmatr, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  double newx = pd->atoms[v].x_orign * curmatr.matr[0][0] + pd->atoms[v].y_orign * curmatr.matr[1][0] + pd->atoms[v].z_orign * curmatr.matr[2][0] + curmatr.matr[3][0];
  double newy = pd->atoms[v].x_orign * curmatr.matr[0][1] + pd->atoms[v].y_orign * curmatr.matr[1][1] + pd->atoms[v].z_orign * curmatr.matr[2][1] + curmatr.matr[3][1];
  double newz = pd->atoms[v].x_orign * curmatr.matr[0][2] + pd->atoms[v].y_orign * curmatr.matr[1][2] + pd->atoms[v].z_orign * curmatr.matr[2][2] + curmatr.matr[3][2];
  pd->atoms[v].x = newx;
  pd->atoms[v].y = newy;
  pd->atoms[v].z = newz;
  for (auto child : tr.children[v]) {
    Matrix4 curmatr_v = curmatr;
    if (!tr.bond_to_parent[child]) {
      rec(pd, tr, child, curmatr_v, encoding_inv, x);
    } else {
      double newx1 = pd->atoms[child].x_orign * curmatr_v.matr[0][0] + pd->atoms[child].y_orign * curmatr_v.matr[1][0] + pd->atoms[child].z_orign * curmatr_v.matr[2][0] + curmatr_v.matr[3][0] - newx;
      double newy1 = pd->atoms[child].x_orign * curmatr_v.matr[0][1] + pd->atoms[child].y_orign * curmatr_v.matr[1][1] + pd->atoms[child].z_orign * curmatr_v.matr[2][1] + curmatr_v.matr[3][1] - newy;
      double newz1 = pd->atoms[child].x_orign * curmatr_v.matr[0][2] + pd->atoms[child].y_orign * curmatr_v.matr[1][2] + pd->atoms[child].z_orign * curmatr_v.matr[2][2] + curmatr_v.matr[3][2] - newz;

      curmatr_v = mlt(curmatr_v, shift_by_vec(-newx, -newy, -newz));
      curmatr_v = mlt(curmatr_v, rotate(newx1, newy1, newz1, x[encoding_inv[child]]));
      curmatr_v = mlt(curmatr_v, shift_by_vec(newx, newy, newz));
      rec(pd, tr, child, curmatr_v, encoding_inv, x);
    }
  }
}

void transformation(hess::Molecule *pd, simplified_tree &tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  Matrix4 init;
  init.matr[0][0] = 1;
  init.matr[1][1] = 1;
  init.matr[2][2] = 1;
  init.matr[3][3] = 1;
  rec(pd, tr, tr.rootid, init, encoding_inv, x);
  Matrix4 M1 = global_rotate_by_x(x[encoding_inv[tr.sz]]);
  Matrix4 M2 = global_rotate_by_y(x[encoding_inv[tr.sz + 1]]);
  Matrix4 M3 = global_rotate_by_z(x[encoding_inv[tr.sz + 2]]);
  Matrix4 M = mlt(M1, M2);
  M = mlt(M, M3);
  for (int i = 0; i < tr.sz; i++) {
    double newx = pd->atoms[i].x * M.matr[0][0] + pd->atoms[i].y * M.matr[1][0] + pd->atoms[i].z * M.matr[2][0] + x[encoding_inv[tr.sz + 3]];
    double newy = pd->atoms[i].x * M.matr[0][1] + pd->atoms[i].y * M.matr[1][1] + pd->atoms[i].z * M.matr[2][1] + x[encoding_inv[tr.sz + 4]];
    double newz = pd->atoms[i].x * M.matr[0][2] + pd->atoms[i].y * M.matr[1][2] + pd->atoms[i].z * M.matr[2][2] + x[encoding_inv[tr.sz + 5]];
    pd->atoms[i].x = newx;
    pd->atoms[i].y = newy;
    pd->atoms[i].z = newz;
  }
}

Matrix4 rotate_by_x_vec_deriv(double angle) {
  Matrix4 res;
  res.matr[0][0] = 0.0;
  res.matr[3][3] = 0.0;
  double sn = sin(angle), cs = cos(angle);
  res.matr[1][1] = -sn, res.matr[1][2] = -cs, res.matr[2][1] = cs, res.matr[2][2] = -sn;
  return res;
}

Matrix4 rotate_deriv(double x, double y, double z, double angle) {
  auto basis = getbasis(x, y, z);
  Matrix4 res = transpose(basis);
  res = mlt(res, rotate_by_x_vec_deriv(angle));
  res = mlt(res, basis);
  return res;
}

void rec_deriv(hess::Molecule *pd, simplified_tree &tr, vector<double>&rotations, int v, Matrix4 curmatr, Matrix4 curder, int alpha_id, int meet_alpha_id, vector<hess::Vec3d>& deriva) {
  double newx = pd->atoms[v].x * curmatr.matr[0][0] + pd->atoms[v].y * curmatr.matr[1][0] + pd->atoms[v].z * curmatr.matr[2][0] + curmatr.matr[3][0];
  double newy = pd->atoms[v].x * curmatr.matr[0][1] + pd->atoms[v].y * curmatr.matr[1][1] + pd->atoms[v].z * curmatr.matr[2][1] + curmatr.matr[3][1];
  double newz = pd->atoms[v].x * curmatr.matr[0][2] + pd->atoms[v].y * curmatr.matr[1][2] + pd->atoms[v].z * curmatr.matr[2][2] + curmatr.matr[3][2];
  double newxd = pd->atoms[v].x * curder.matr[0][0] + pd->atoms[v].y * curder.matr[1][0] + pd->atoms[v].z * curder.matr[2][0] + curder.matr[3][0];
  double newyd = pd->atoms[v].x * curder.matr[0][1] + pd->atoms[v].y * curder.matr[1][1] + pd->atoms[v].z * curder.matr[2][1] + curder.matr[3][1];
  double newzd = pd->atoms[v].x * curder.matr[0][2] + pd->atoms[v].y * curder.matr[1][2] + pd->atoms[v].z * curder.matr[2][2] + curder.matr[3][2];
  if (meet_alpha_id) {
    deriva[v].x = newxd;
    deriva[v].y = newyd;
    deriva[v].z = newzd;
  } else {
    deriva[v].x = 0;
    deriva[v].y = 0;
    deriva[v].z = 0;
    pd->atoms[v].x = 0;
    pd->atoms[v].y = 0;
    pd->atoms[v].z = 0;
  }
  if (v == alpha_id) {
    meet_alpha_id = 1;

  }

  for (auto child : tr.children[v]) {
    if (!tr.bond_to_parent[child]) {
      rec_deriv(pd, tr, rotations, child, curmatr, curder, alpha_id, meet_alpha_id, deriva);
    } else {
      double newx1 = pd->atoms[child].x * curmatr.matr[0][0] + pd->atoms[child].y * curmatr.matr[1][0] + pd->atoms[child].z * curmatr.matr[2][0] + curmatr.matr[3][0] - newx;
      double newy1 = pd->atoms[child].x * curmatr.matr[0][1] + pd->atoms[child].y * curmatr.matr[1][1] + pd->atoms[child].z * curmatr.matr[2][1] + curmatr.matr[3][1] - newy;
      double newz1 = pd->atoms[child].x * curmatr.matr[0][2] + pd->atoms[child].y * curmatr.matr[1][2] + pd->atoms[child].z * curmatr.matr[2][2] + curmatr.matr[3][2] - newz;
      curmatr = mlt(curmatr, shift_by_vec(-newx, -newy, -newz));
      curmatr = mlt(curmatr, rotate(newx1, newy1, newz1, rotations[child]));
      curmatr = mlt(curmatr, shift_by_vec(newx, newy, newz));
      double newx2 = pd->atoms[child].x * curder.matr[0][0] + pd->atoms[child].y * curder.matr[1][0] + pd->atoms[child].z * curder.matr[2][0] + curder.matr[3][0] - newx;
      double newy2 = pd->atoms[child].x * curder.matr[0][1] + pd->atoms[child].y * curder.matr[1][1] + pd->atoms[child].z * curder.matr[2][1] + curder.matr[3][1] - newy;
      double newz2 = pd->atoms[child].x * curder.matr[0][2] + pd->atoms[child].y * curder.matr[1][2] + pd->atoms[child].z * curder.matr[2][2] + curder.matr[3][2] - newz;
      curder = mlt(curder, shift_by_vec(-newx, -newy, -newz));
      if ((child == alpha_id) && tr.bond_to_parent[alpha_id]) {
        curder = mlt(curder, rotate_deriv(newx2, newy2, newz2, rotations[child]));
      } else curder = mlt(curder, rotate(newx2, newy2, newz2, rotations[child]));
      curder = mlt(curder, shift_by_vec(newx, newy, newz));
      rec_deriv(pd, tr, rotations, child, curmatr, curder, alpha_id, meet_alpha_id, deriva);
    }
  }
}


