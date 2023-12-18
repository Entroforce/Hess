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

#include "fast_gradient.h"
#include "optimization_methods.h"

double fast_gradient(ScoringTerms& st, Optimizable_molecule* mol, const Eigen::VectorXd& x, Eigen::VectorXd& grad, const vector<int>& encoding_inv) {
  int lig_atoms = mol->ligand->get_atoms_count();
  int rots_count = mol->rot_bonds_count;
  double max_cutoff_sqr = st.get_max_cutoff_sqr();
  const int n = num_atom_types();
  thread_local static vector<hess::Vec3d> da_dalpha;
  if (!mol->type) {
    assert(rots_count + 6 == grad.size());
    da_dalpha.resize((rots_count + 3) * lig_atoms);
    calc_all_da_dalpha(da_dalpha, mol->ligand, mol->tr, encoding_inv, x, mol->rotsmap); //transformation inside
  } else if (mol->type == 1) {
    assert(6 == grad.size());
    da_dalpha.resize(3 * lig_atoms);
    calc_all_da_dalpha_type1(da_dalpha, mol->ligand, mol->tr, encoding_inv, x);
  }
  double energy = 0.0;
  double e_intra = 0.0;
  for (int a_id = mol->ligand->vertexBegin(); a_id != mol->ligand->vertexEnd(); a_id = mol->ligand->vertexNext(a_id)) {
    hess::Atom* a = mol->ligand->get_atom(a_id);
    smt t1 = a->atom_type;
    if (t1 >= n || is_hydrogen(t1)) {
      continue;
    }
    const Array3D<double> &mgn = mol->grid_new[t1];
    hess::Vec3d buf_diff;
    hess::Vec3d prod;
    hess::Vec3d location(a->get_vector());
    buf_diff.diff(location, mol->init_vec);
    prod.elementwise_product(buf_diff, mol->factor_vec);
    int region[3] = {0};
    hess::Vec3d miss;
    int a_0[3] = {0};
    for (int i = 0; i < 3; i++) {
      if (prod[i] < 0) {
        miss[i] = -prod[i];
        region[i] = -1;
        a_0[i] = 0;
        prod[i] = 0;
      } else if (prod[i] >= mol->dim_fl_minus_1[i]) {
        miss[i] = prod[i] - mol->dim_fl_minus_1[i];
        region[i] = 1;
        assert(mgn.dim(i) >= 2);
        a_0[i] = mgn.dim(i) - 2;
        prod[i] = 1;
      } else {
        region[i] = 0;
        a_0[i] = std::size_t(prod[i]);
        prod[i] = prod[i] - a_0[i];
      }
      if (prod[i] < 0 || prod[i] > 1 || a_0[i] < 0 || a_0[i] + 1 >= mgn.dim(i)) {
        fprintf(stderr, "\na_0 = %lf, prod = %lf\n", a_0[i], prod[i]);
      }
      assert(prod[i] >= 0);
      assert(prod[i] <= 1);
      assert(a_0[i] >= 0);
      assert(a_0[i] + 1 < mgn.dim(i));
    }
    const double penalty = mol->slope * (hess::Vec3d::dot(miss, mol->factor_inv));
    assert(penalty > -epsilon_fl);
    const int x0 = a_0[0];
    const int y0 = a_0[1];
    const int z0 = a_0[2];
    const int x1 = x0 + 1;
    const int y1 = y0 + 1;
    const int z1 = z0 + 1;
    const double f000 = mgn(x0, y0, z0);
    const double f100 = mgn(x1, y0, z0);
    const double f010 = mgn(x0, y1, z0);
    const double f110 = mgn(x1, y1, z0);
    const double f001 = mgn(x0, y0, z1);
    const double f101 = mgn(x1, y0, z1);
    const double f011 = mgn(x0, y1, z1);
    const double f111 = mgn(x1, y1, z1);
    const double x = prod.x;
    const double y = prod.y;
    const double z = prod.z;
    const double mx = 1 - x;
    const double my = 1 - y;
    const double mz = 1 - z;
    double f =
            f000 * mx * my * mz +
            f100 * x * my * mz +
            f010 * mx * y * mz +
            f110 * x * y * mz +
            f001 * mx * my * z +
            f101 * x * my * z +
            f011 * mx * y * z +
            f111 * x * y * z;
    double x_g = 0.0, y_g = 0.0, z_g = 0.0;
    if (mol->grid_with_derivs) {
      Array3D<hess::Vec3d>& mgdn = mol->grid_deriv_new[t1];
      const hess::Vec3d f000_d = mgdn(x0, y0, z0);
      const hess::Vec3d f100_d = mgdn(x1, y0, z0);
      const hess::Vec3d f010_d = mgdn(x0, y1, z0);
      const hess::Vec3d f110_d = mgdn(x1, y1, z0);
      const hess::Vec3d f001_d = mgdn(x0, y0, z1);
      const hess::Vec3d f101_d = mgdn(x1, y0, z1);
      const hess::Vec3d f011_d = mgdn(x0, y1, z1);
      const hess::Vec3d f111_d = mgdn(x1, y1, z1);
      x_g =
              f000_d.x * mx * my * mz +
              f100_d.x * x * my * mz +
              f010_d.x * mx * y * mz +
              f110_d.x * x * y * mz +
              f001_d.x * mx * my * z +
              f101_d.x * x * my * z +
              f011_d.x * mx * y * z +
              f111_d.x * x * y * z;

      y_g =
              f000_d.y * mx * my * mz +
              f100_d.y * x * my * mz +
              f010_d.y * mx * y * mz +
              f110_d.y * x * y * mz +
              f001_d.y * mx * my * z +
              f101_d.y * x * my * z +
              f011_d.y * mx * y * z +
              f111_d.y * x * y * z;

      z_g =
              f000_d.z * mx * my * mz +
              f100_d.z * x * my * mz +
              f010_d.z * mx * y * mz +
              f110_d.z * x * y * mz +
              f001_d.z * mx * my * z +
              f101_d.z * x * my * z +
              f011_d.z * mx * y * z +
              f111_d.z * x * y * z;
    }
    else {
      x_g =
              f000 * (-1) * my * mz +
              f100 * 1 * my * mz +
              f010 * (-1) * y * mz +
              f110 * 1 * y * mz +
              f001 * (-1) * my * z +
              f101 * 1 * my * z +
              f011 * (-1) * y * z +
              f111 * 1 * y * z;

      y_g =
              f000 * mx * (-1) * mz +
              f100 * x * (-1) * mz +
              f010 * mx * 1 * mz +
              f110 * x * 1 * mz +
              f001 * mx * (-1) * z +
              f101 * x * (-1) * z +
              f011 * mx * 1 * z +
              f111 * x * 1 * z;

      z_g =
              f000 * mx * my * (-1) +
              f100 * x * my * (-1) +
              f010 * mx * y * (-1) +
              f110 * x * y * (-1) +
              f001 * mx * my * 1 +
              f101 * x * my * 1 +
              f011 * mx * y * 1 +
              f111 * x * y * 1;
    }
    hess::Vec3d gradient(x_g * mol->factor_vec[0], y_g * mol->factor_vec[1], z_g * mol->factor_vec[2]);
    curl(f, gradient, mol->cap[0]);
    hess::Vec3d deriv_a_dtdr((region[0] == 0) ? gradient.x : mol->slope * region[0],
            (region[1] == 0) ? gradient.y : mol->slope * region[1],
            (region[2] == 0) ? gradient.z : mol->slope * region[2]);
    energy += (f + penalty);
    if (!mol->type) {
      double dot_dr_da_ori_x = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[mol->tr.sz * rots_count + a_id]);
      double dot_dr_da_ori_y = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[mol->tr.sz * (rots_count + 1) + a_id]);
      double dot_dr_da_ori_z = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[mol->tr.sz * (rots_count + 2) + a_id]);
      grad[rots_count] += dot_dr_da_ori_x;
      grad[rots_count + 1] += dot_dr_da_ori_y;
      grad[rots_count + 2] += dot_dr_da_ori_z;
      grad[rots_count + 3] += deriv_a_dtdr.x;
      grad[rots_count + 4] += deriv_a_dtdr.y;
      grad[rots_count + 5] += deriv_a_dtdr.z;
      for (int tor_id = 0; tor_id < rots_count; tor_id++) {
        const hess::Vec3d& deriv_a = da_dalpha[mol->tr.sz * tor_id + a_id];
        double dot_dr_da = hess::Vec3d::dot(deriv_a_dtdr, deriv_a);
        grad[tor_id] += dot_dr_da;
      }
    } else if (mol->type == 1) {
      double dot_dr_da_ori_x = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[a_id]);
      double dot_dr_da_ori_y = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[mol->tr.sz + a_id]);
      double dot_dr_da_ori_z = hess::Vec3d::dot(deriv_a_dtdr, da_dalpha[mol->tr.sz * 2 + a_id]);
      grad[0] += dot_dr_da_ori_x;
      grad[1] += dot_dr_da_ori_y;
      grad[2] += dot_dr_da_ori_z;
      grad[3] += deriv_a_dtdr.x;
      grad[4] += deriv_a_dtdr.y;
      grad[5] += deriv_a_dtdr.z;
    }
  }
  //Start calc intramolecular
  if (!mol->type) {
    thread_local static vector<hess::Vec3d> pairs_precal_plus, pairs_precal_minus;
    pairs_precal_plus.resize(mol->ligand->get_atoms_count());
    pairs_precal_minus.resize(mol->ligand->get_atoms_count());
    memset(&pairs_precal_plus[0], 0, sizeof (hess::Vec3d) * pairs_precal_plus.size());
    memset(&pairs_precal_minus[0], 0, sizeof (hess::Vec3d) * pairs_precal_minus.size());
    for (auto &pair : mol->ligand->pairs) {
      int a_id = pair.first;
      int b_id = pair.second;
      hess::Atom* a = mol->ligand->get_atom(a_id);
      hess::Atom* b = mol->ligand->get_atom(b_id);
      smt t1 = a->atom_type;
      smt t2 = b->atom_type;
      hess::Vec3d a_coord(a->x, a->y, a->z);
      hess::Vec3d b_coord(b->x, b->y, b->z);
      hess::Vec3d diff_buf;
      diff_buf.diff(a_coord, b_coord);
      double r2 = diff_buf.lengthSqr();
      if (r2 < max_cutoff_sqr) {
        double r = sqrt(r2);
        double this_deriv[3] = {0};
        st.calc_term_deriv_intramolecular_without_a(this_deriv, t1, t2, r, &da_dalpha[0],
                mol->tr.sz * (rots_count) + a_id, mol->tr.sz * (rots_count + 1) + a_id, mol->tr.sz * (rots_count + 2) + a_id,
                mol->tr.sz * (rots_count) + b_id, mol->tr.sz * (rots_count + 1) + b_id, mol->tr.sz * (rots_count + 2) + b_id,
                diff_buf);
        grad[rots_count] += this_deriv[0];
        grad[rots_count + 1] += this_deriv[1];
        grad[rots_count + 2] += this_deriv[2];
        e_intra += st.eval_fast(t1, t2, r);
        hess::Vec3d d_dterms_mul = st.calc_term_deriv_part_dist(a_coord, b_coord);
        double d_terms = st.calc_term_deriv_part_terms(t1, t2, r);
        d_dterms_mul.x *= d_terms;
        d_dterms_mul.y *= d_terms;
        d_dterms_mul.z *= d_terms;
        pairs_precal_plus[a_id].add(d_dterms_mul);
        pairs_precal_minus[b_id].add(d_dterms_mul);
      }
    }
    for (int tor_id = 0; tor_id < rots_count; tor_id++) {
      double acc_torsion_deriv = 0.0;
      for (int a_id = 0; a_id < mol->ligand->get_atoms_count(); a_id++) {
        const hess::Vec3d &deriv_a = da_dalpha[mol->tr.sz * tor_id + a_id];
        acc_torsion_deriv += hess::Vec3d::dot(pairs_precal_plus[a_id], deriv_a);
        acc_torsion_deriv -= hess::Vec3d::dot(pairs_precal_minus[a_id], deriv_a);
      }
      grad[tor_id] += acc_torsion_deriv;
    }
  }    //1e-13 max diff
  else if (mol->type == 1) {
    for (auto &pair : mol->ligand->pairs) {
      int a_id = pair.first;
      int b_id = pair.second;
      hess::Atom* a = mol->ligand->get_atom(a_id);
      hess::Atom* b = mol->ligand->get_atom(b_id);
      smt t1 = a->atom_type;
      smt t2 = b->atom_type;
      hess::Vec3d diff_buf(a->x - b->x, a->y - b->y, a->z - b->z);
      double r2 = diff_buf.lengthSqr();
      if (r2 < max_cutoff_sqr) {
        double r = sqrt(r2);
        double this_deriv[3] = {0};
        st.calc_term_deriv_intramolecular_type_1(this_deriv, t1, t2, r, &da_dalpha[0],
                a_id, mol->tr.sz + a_id, mol->tr.sz * 2 + a_id,
                b_id, mol->tr.sz + b_id, mol->tr.sz * 2 + b_id, diff_buf);
        grad[0] += this_deriv[0];
        grad[1] += this_deriv[1];
        grad[2] += this_deriv[2];
        e_intra += st.eval_fast(t1, t2, r);

      }
    }
  }
  energy += e_intra;
  //END calc intramolecular
  return energy;
}