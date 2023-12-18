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

#include "features.h"

class ScoringTerms {
public:
  virtual double get_max_cutoff_sqr();
  virtual double eval_fast(smt t1, smt t2, double r);
  virtual double conf_independent(hess::Molecule* ligand, double e) = 0;
  virtual double conf_independent(hess::Molecule* ligand, ConfIndependentInputs& in, double e) = 0;
  virtual void calc_term_deriv_intramolecular_without_a(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) = 0;
  virtual double calc_term_deriv_intramolecular_deriv_a_only(smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d& diff) = 0;
  virtual void calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) = 0;
  virtual void calc_term_deriv_intramolecular_type_1(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff) = 0;
  virtual void calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d* deriv_ori_a, hess::Vec3d* deriv_ori_b, hess::Vec3d& diff) = 0;
  virtual hess::Vec3d calc_term_deriv_part_dist(hess::Vec3d& v1, hess::Vec3d& v2) = 0;
  virtual double calc_term_deriv_part_terms(smt t1, smt t2, double r) = 0;
  virtual double get_weight(int id);
  virtual int size();
  virtual ~ScoringTerms();
};

class VinardoWeightedTerms : public ScoringTerms {
public:

  VinardoWeightedTerms();
  ~VinardoWeightedTerms();
  int size();
  double get_weight(int id);
  double calc_max_cutoff();
  double get_max_cutoff_sqr();
  double get_max_cutoff();
  double calc_term_deriv_intramolecular_deriv_a_only(smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b, hess::Vec3d& diff);
  void calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b,
          hess::Vec3d *deriv, int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff);
  void calc_term_deriv_intramolecular(double* grad, smt t1, smt t2, double r, hess::Vec3d& deriv_a, hess::Vec3d& deriv_b,
          hess::Vec3d* deriv_ori_a, hess::Vec3d* deriv_ori_b, hess::Vec3d& diff);
  void calc_term_deriv_intramolecular_without_a(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv,
          int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff);
  void calc_term_deriv_intramolecular_type_1(double* grad, smt t1, smt t2, double r, hess::Vec3d *deriv,
          int ori_a0, int ori_a1, int ori_a2, int ori_b0, int ori_b1, int ori_b2, hess::Vec3d& diff);
  hess::Vec3d calc_term_deriv_part_dist(hess::Vec3d& v1, hess::Vec3d& v2);
  double calc_term_deriv_part_terms(smt t1, smt t2, double r);
  double eval_fast(smt t1, smt t2, double r);
  double conf_independent(hess::Molecule* ligand, ConfIndependentInputs& in, double e);
  double conf_independent(hess::Molecule* ligand, double e);

private:
  unsigned terms_count;
  double max_cutoff;
  double max_cutoff_sqr;
  Gauss g;
  Repulsion r;
  Hydrophobic h;
  Non_dir_h_bond nd;
  Num_tors_div nt;
  double weights[5] = {-0.045, 0.8, -0.035, -0.6, 5 * 0.02 / 0.1 - 1};
};


extern VinardoWeightedTerms vinardo_sc;

double calc_intramolecular_energy(hess::Molecule *lig, hess::Molecule *rec, const char *scoring);
double eval_adjusted(ScoringTerms &st, hess::Molecule *ligand, hess::Molecule *receptor);
double calc_affinity(hess::Molecule* lig, hess::Molecule* rec, const char *scoring);
double calc_affinity_with_confind(hess::Molecule *lig, hess::Molecule *rec, ConfIndependentInputs& in, const char *scoring);
bool geometry_changed(hess::Molecule* mol);

inline bool not_max(double x) {
  return (x < 0.1 * max_fl);
}

inline void curl(double &e, double v = 1000) {
  if (e > 0 && not_max(v)) {
    double tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
    e *= tmp;
  }
}

inline void curl(double &e, hess::Vec3d& deriv, double v = 1000) {
  if (e > 0 && not_max(v)) {
    double tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
    e *= tmp;
    deriv.scale(sqr(tmp));
  }
}