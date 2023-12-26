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

/**********************************************************************
 * Parts of this file were borrowed from smina repository
 * published under Apache License 2.0
 * https://github.com/mwojcikowski/smina/blob/master/src/lib/model.cpp
 * Copyright (c) 2006-2010, The Scripps Research Institute 
 ***********************************************************************/

#include "scoring.h"

using namespace hess;

double ScoringTerms::get_max_cutoff_sqr() {
  return 0.0;
}

double ScoringTerms::eval_fast(smt t1, smt t2, double r) {
  return 0.0;
}

double ScoringTerms::get_weight(int id) {
  return 0.0;
}

int ScoringTerms::size() {
  return 0;
}

ScoringTerms::~ScoringTerms() {
}

VinardoWeightedTerms::VinardoWeightedTerms() {
  g = Gauss(0, 0.8, 8);
  r = Repulsion(0, 8);
  h = Hydrophobic(0.0, 2.5, 8);
  nd = Non_dir_h_bond(-0.6, 0, 8);
  nt = Num_tors_div();
  terms_count = 5;
  this->max_cutoff = 8;
  this->max_cutoff_sqr = 64;
}

VinardoWeightedTerms::~VinardoWeightedTerms() {
}

int VinardoWeightedTerms::size() {
  return this->terms_count;
}

double VinardoWeightedTerms::get_weight(int id) {
  assert(id < terms_count);
  return this->weights[id];
}

double VinardoWeightedTerms::get_max_cutoff_sqr() {
  return this->max_cutoff_sqr;
}

double VinardoWeightedTerms::get_max_cutoff() {
  return this->max_cutoff;
}

double VinardoWeightedTerms::calc_term_deriv_intramolecular_deriv_a_only(smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d& diff) {
  //        if (!diff.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  diff.normalize();
  hess::Vec3d diff_deriv;
  diff_deriv.diff(deriv_a, deriv_b);
  double deriv_res = hess::Vec3d::dot(diff, diff_deriv);
  double dterm_acc = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_res +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_res +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_res +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_res;
  return dterm_acc;
}

void VinardoWeightedTerms::calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) {
  //        if (!diff.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  diff.normalize();
  hess::Vec3d diff_deriv;
  hess::Vec3d diff_ori_x;
  hess::Vec3d diff_ori_y;
  hess::Vec3d diff_ori_z;
  diff_deriv.diff(deriv_a, deriv_b);
  diff_ori_x.diff(deriv[ori_a0], deriv[ori_b0]);
  diff_ori_y.diff(deriv[ori_a1], deriv[ori_b1]);
  diff_ori_z.diff(deriv[ori_a2], deriv[ori_b2]);
  double deriv_ori_x = hess::Vec3d::dot(diff, diff_ori_x);
  double deriv_ori_y = hess::Vec3d::dot(diff, diff_ori_y);
  double deriv_ori_z = hess::Vec3d::dot(diff, diff_ori_z);
  double deriv_res = hess::Vec3d::dot(diff, diff_deriv);

  double dterm_acc = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_res +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_res +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_res +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_res;

  double dori_acc_x = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_x +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_x +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_x +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_x;

  double dori_acc_y = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_y +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_y +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_y +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_y;

  double dori_acc_z = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_z +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_z +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_z +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_z;
  grad[0] = dterm_acc;
  grad[1] = dori_acc_x;
  grad[2] = dori_acc_y;
  grad[3] = dori_acc_z;
}

void VinardoWeightedTerms::calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d* deriv_ori_a, hess::Vec3d* deriv_ori_b, hess::Vec3d& diff) {
  //        if (!diff.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  diff.normalize();
  hess::Vec3d diff_deriv;
  hess::Vec3d diff_ori_x;
  hess::Vec3d diff_ori_y;
  hess::Vec3d diff_ori_z;
  diff_deriv.diff(deriv_a, deriv_b);
  diff_ori_x.diff(deriv_ori_a[0], deriv_ori_b[0]);
  diff_ori_y.diff(deriv_ori_a[1], deriv_ori_b[1]);
  diff_ori_z.diff(deriv_ori_a[2], deriv_ori_b[2]);
  double deriv_ori_x = hess::Vec3d::dot(diff, diff_ori_x);
  double deriv_ori_y = hess::Vec3d::dot(diff, diff_ori_y);
  double deriv_ori_z = hess::Vec3d::dot(diff, diff_ori_z);
  double deriv_res = hess::Vec3d::dot(diff, diff_deriv);

  double dterm_acc = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_res +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_res +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_res +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_res;

  double dori_acc_x = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_x +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_x +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_x +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_x;

  double dori_acc_y = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_y +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_y +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_y +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_y;

  double dori_acc_z = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_z +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_z +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_z +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_z;
  grad[0] = dterm_acc;
  grad[1] = dori_acc_x;
  grad[2] = dori_acc_y;
  grad[3] = dori_acc_z;
}

void VinardoWeightedTerms::calc_term_deriv_intramolecular_without_a(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) {
  //        if (!diff.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  diff.normalize();
  hess::Vec3d diff_ori_x;
  hess::Vec3d diff_ori_y;
  hess::Vec3d diff_ori_z;
  diff_ori_x.diff(deriv[ori_a0], deriv[ori_b0]);
  diff_ori_y.diff(deriv[ori_a1], deriv[ori_b1]);
  diff_ori_z.diff(deriv[ori_a2], deriv[ori_b2]);
  double deriv_ori_x = hess::Vec3d::dot(diff, diff_ori_x);
  double deriv_ori_y = hess::Vec3d::dot(diff, diff_ori_y);
  double deriv_ori_z = hess::Vec3d::dot(diff, diff_ori_z);

  double dori_acc_x = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_x +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_x +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_x +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_x;

  double dori_acc_y = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_y +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_y +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_y +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_y;

  double dori_acc_z = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_z +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_z +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_z +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_z;
  grad[0] = dori_acc_x;
  grad[1] = dori_acc_y;
  grad[2] = dori_acc_z;
}

void VinardoWeightedTerms::calc_term_deriv_intramolecular_type_1(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) {
  //        if (!diff.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  diff.normalize();
  hess::Vec3d diff_ori_x;
  hess::Vec3d diff_ori_y;
  hess::Vec3d diff_ori_z;
  diff_ori_x.diff(deriv[ori_a0], deriv[ori_b0]);
  diff_ori_y.diff(deriv[ori_a1], deriv[ori_b1]);
  diff_ori_z.diff(deriv[ori_a2], deriv[ori_b2]);
  double deriv_ori_x = hess::Vec3d::dot(diff, diff_ori_x);
  double deriv_ori_y = hess::Vec3d::dot(diff, diff_ori_y);
  double deriv_ori_z = hess::Vec3d::dot(diff, diff_ori_z);

  double dori_acc_x = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_x +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_x +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_x +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_x;

  double dori_acc_y = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_y +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_y +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_y +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_y;

  double dori_acc_z = this->g.eval_deriv(t1, t2, r) * this->weights[0] * deriv_ori_z +
          this->r.eval_deriv(t1, t2, r) * this->weights[1] * deriv_ori_z +
          this->h.eval_deriv(t1, t2, r) * this->weights[2] * deriv_ori_z +
          this->nd.eval_deriv(t1, t2, r) * this->weights[3] * deriv_ori_z;
  grad[0] = dori_acc_x;
  grad[1] = dori_acc_y;
  grad[2] = dori_acc_z;
}

hess::Vec3d VinardoWeightedTerms::calc_term_deriv_part_dist(hess::Vec3d& v1, hess::Vec3d& v2) {
  hess::Vec3d rd;
  rd.diff(v1, v2);
  //        if (!rd.normalize())
  //            throw HessException("Normalization error for derivative calculation.\n");
  rd.normalize();
  return rd;
}

double VinardoWeightedTerms::calc_term_deriv_part_terms(smt t1, smt t2, double r) {
  double dt_acc =
          this->g.eval_deriv(t1, t2, r) * this->weights[0]
          + this->r.eval_deriv(t1, t2, r) * this->weights[1]
          + this->h.eval_deriv(t1, t2, r) * this->weights[2]
          + this->nd.eval_deriv(t1, t2, r) * this->weights[3]
          ;
  return dt_acc;
}

double VinardoWeightedTerms::eval_fast(smt t1, smt t2, double r) { // intentionally not checking for cutoff
  double acc =
          this->g.eval(t1, t2, r) * this->weights[0]
          + this->r.eval(t1, t2, r) * this->weights[1]
          + this->h.eval(t1, t2, r) * this->weights[2]
          + this->nd.eval(t1, t2, r) * this->weights[3]
          ;
  return acc;

}

double VinardoWeightedTerms::conf_independent(hess::Molecule* ligand, ConfIndependentInputs& in, double e) {
  Term &t = nt;
  double weight = this->get_weight(4);
  e = t.eval(in, e, weight);
  return e;
}

double VinardoWeightedTerms::conf_independent(hess::Molecule* ligand, double e) {
  ConfIndependentInputs in(ligand);
  Term &t = nt;
  double weight = this->get_weight(4);
  e = t.eval(in, e, weight);
  return e;
}

VinardoWeightedTerms vinardo_sc = VinardoWeightedTerms();

double eval_intramolecular_energy(ScoringTerms &st, hess::Molecule *ligand, hess::Molecule *receptor) {
  double e_intra = 0.0;
  double max_cutoff_sqr = st.get_max_cutoff_sqr();
  const vector<pair<int, int>> &lig_pairs = ligand->get_lig_pairs();
  int i = 0;
  for (auto pair : lig_pairs) {
    hess::Atom* a = ligand->get_atom(pair.first);
    hess::Atom* b = ligand->get_atom(pair.second);
    double r2 = distance_sqr(*a, *b);
    if (r2 < max_cutoff_sqr) {
      double tmp = st.eval_fast(a->atom_type, b->atom_type, sqrt(r2));
      curl(tmp);
      e_intra += tmp;
//      printf("Pair %d: %d, %d, %.4g, %.4g\n", i++, pair.first+1, pair.second+1, distance(*a, *b), tmp);
    }
  }
  return e_intra;
}

double eval_adjusted(hess::Molecule *ligand, hess::Molecule *receptor) {
  double e = 0.0;
  double max_cutoff_sqr = vinardo_sc.get_max_cutoff_sqr();
  const int n = num_atom_types();
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    Atom *a = ligand->get_atom(a_id);
    double this_e = 0.0;
    smt t1 = a->atom_type;
    if (t1 >= n || is_hydrogen(t1)) {
      continue;
    }
    for (int b_id = receptor->vertexBegin(); b_id != receptor->vertexEnd(); b_id = receptor->vertexNext(b_id)) {
      Atom *b = receptor->get_atom(b_id);
      smt t2 = b->atom_type;
      if (t2 >= n || is_hydrogen(t2)) {
        continue;
      }
      double r2 = distance_sqr(*a, *b);
      if (r2 < max_cutoff_sqr) {
        this_e += vinardo_sc.eval_fast(t1, t2, sqrt(r2));
      }
    }
    curl(this_e);
    e += this_e;
  }
  double e_conf_independet = vinardo_sc.conf_independent(ligand, e);
  return e_conf_independet;
}

double eval_adjusted_upd_with_confind(ScoringTerms &st, hess::Molecule *ligand, hess::Molecule *receptor, ConfIndependentInputs& in) {
  double e = 0.0;
  double max_cutoff_sqr = st.get_max_cutoff_sqr();
  const int n = num_atom_types();
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    Atom *a = ligand->get_atom(a_id);
    double this_e = 0.0;
    smt t1 = a->atom_type;
    if (t1 >= n || is_hydrogen(t1)) {
      continue;
    }
    for (int b_id = receptor->vertexBegin(); b_id != receptor->vertexEnd(); b_id = receptor->vertexNext(b_id)) {
      Atom *b = receptor->get_atom(b_id);
      smt t2 = b->atom_type;
      if (t2 >= n || is_hydrogen(t2)) {
        continue;
      }
      double r2 = distance_sqr(*a, *b);
      if (r2 < max_cutoff_sqr) {
        this_e += st.eval_fast(a->atom_type, b->atom_type, sqrt(r2));
      }
    }
    curl(this_e);
    e += this_e;
  }
  return st.conf_independent(ligand, in, e);
}

bool geometry_changed(hess::Molecule* mol) {
  set <pair<int, int>> ex_bonds;
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    for (int b_id = mol->vertexBegin(); b_id != mol->vertexEnd(); b_id = mol->vertexNext(b_id)) {
      hess::Atom *a = mol->get_atom(a_id);
      hess::Atom *b = mol->get_atom(b_id);
      if ((a_id == b_id) || (ex_bonds.find(make_pair(b_id, a_id)) != ex_bonds.end())) {
        continue;
      }
      if (mol->getAtomNumber(a_id) == 1 || mol->getAtomNumber(b_id) == 1) {
        continue;
      }

      int bond_order = mol->determine_order(a_id, b_id);
      if (bond_order) {
        int old_bond_id = mol->findEdgeIndex(a_id, b_id);
        if (old_bond_id == -1)
          return true;
      } else {
        int old_bond_id = mol->findEdgeIndex(a_id, b_id);
        if (old_bond_id != -1)
          return true;
      }
      ex_bonds.insert(make_pair(a_id, b_id));
    }
  }
  return false;
}

double calc_intramolecular_energy(hess::Molecule *lig, hess::Molecule *rec) {
  double e_s = 0.0;
  e_s = eval_intramolecular_energy(vinardo_sc, lig, rec);
  return e_s;
}

double calc_affinity_with_confind(hess::Molecule *lig, hess::Molecule *rec, ConfIndependentInputs& in, const char *scoring) {
  double e_s = 0.0;
  e_s = eval_adjusted_upd_with_confind(vinardo_sc, lig, rec, in);
  return e_s;
}

double calc_affinity(hess::Molecule *lig, hess::Molecule *rec) {
  double e_s = 0.0;
  e_s = eval_adjusted(lig, rec);
  return e_s;
}
