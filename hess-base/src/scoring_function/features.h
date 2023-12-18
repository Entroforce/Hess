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

#include "bond_order_determination.h"
#include "model/molecule.h"
#include "math/algebra.h"


void assign_atom_types(hess::Molecule *mol, const string& scoring_type_str);
bool atom_is_Hbond_acceptor(hess::Molecule *mol, int atom_id);
bool bond_is_ester(hess::Molecule *mol, int bgn_id, int end_id);
bool bond_is_carbonyl(hess::Molecule *mol, int bgn_id, int end_id);
bool atom_is_sulfone_oxygen(hess::Molecule *mol, int atom_id);
bool atom_is_nitro_oxygen(hess::Molecule *mol, int atom_id);
bool bond_is_amide(hess::Molecule *mol, int bond_id);
bool is_suppressible_hydrogen(hess::Molecule* mol, int atom_id);
bool is_rot_bond(hess::Molecule *mol, int bond_id);


double gaussian(double x, double width);
double smooth_div(double x, double y);
double slope_step(double x_bad, double x_good, double x);
double slope_step_deriv(double x_bad, double x_good, double x);
double optimal_distance(smt xs_t1, smt xs_t2);
bool bonded_to_heteroatom(hess::Molecule *mol, const int atom_id);
bool bonded_to_HD(hess::Molecule *mol, const int atom_id);
smt adjust_atom_type(smt t, bool Hbonded, bool heteroBonded);
void adjust_atom_types(hess::Molecule *mol);

enum TermType {
  BaseTerm, ConfIndependent, InterMolecular, Additive, DistanceAdditive, ChargeDependent, ChargeIndependent, LastTermType
};

class ConfIndependentInputs {
public:
  ConfIndependentInputs(hess::Molecule* ligand);
  ConfIndependentInputs();
  double num_tors;
  double num_rotors;
  double num_heavy_atoms;
  double num_hydrophobic_atoms;
  double ligand_max_num_h_bonds;
  double num_ligands;
  double ligand_lengths_sum;

private:
  unsigned num_bonded_heavy_atoms(hess::Molecule* ligand, int a_id);
  unsigned atom_rotors(hess::Molecule* ligand, int a_id);
};

class Term {
public:

  Term(double cut = 0.0) : cutoff(cut) {
  }

  virtual ~Term() {
  }

  virtual TermType get_term_type() {
    return TermType::BaseTerm;
  }

  double get_cutoff() {
    return cutoff;
  }

  virtual double eval(smt t1, smt t2, double r) {
    return 0.0;
  }

  virtual double eval_deriv(const smt t1, const smt t2, double r) {
    return 0.0;
  }

  virtual double eval(ConfIndependentInputs &in, double x, double w) {
    return 0.0;
  }

  virtual Term * copy() = 0;
  string term_name;
  double cutoff;
};

class Gauss : public Term {
public:
  TermType type;
  double offset;
  double width;

  Gauss(double offset_ = 0, double width_ = 0.5, double cutoff_ = 8) :
  Term(cutoff_), offset(offset_), width(width_), type(TermType::ChargeIndependent) {
    term_name = string("gauss(o=") + to_string(offset) + ",_w="
            + to_string(width) + ",_c=" + to_string(cutoff) + ")";
  }

  ~Gauss() {
  }

  TermType get_term_type() {
    return type;
  }

  double eval_deriv(const smt t1, const smt t2, double r) {
    double opt = optimal_distance(t1, t2);
    double der = -2 * exp(-sqr((r - opt + offset) / width)) * ((r - opt + offset) / width) / width;
    return der;
  }

  double eval(smt t1, smt t2, double r) {
    return gaussian(r - (optimal_distance(t1, t2) + offset), width);
  }

  Term *copy() {
    return new Gauss(*this);
  }

};

class Repulsion : public Term {
public:
  double offset;
  TermType type;

  Repulsion(double offset_ = 0, double cutoff_ = 8) :
  Term(cutoff_), offset(offset_), type(TermType::ChargeIndependent) {
    term_name = string("repulsion(o=") + to_string(offset) + ",_c="
            + to_string(cutoff) + ")";
  }

  ~Repulsion() {
  }

  TermType get_term_type() {
    return type;
  }

  double eval(smt t1, smt t2, double r) {
    double d = r - (optimal_distance(t1, t2) + offset);
    if (d > 0)
      return 0;
    return d * d;
  }

  double eval_deriv(const smt t1, const smt t2, double r) {
    double d = r - (optimal_distance(t1, t2) + offset);
    if (d > 0)
      return 0;
    return 2 * d;
  }

  Term *copy() {
    return new Repulsion(*this);
  }
};

class Hydrophobic : public Term {
public:
  TermType type;
  double good;
  double bad;

  Hydrophobic(double good_ = 0.5, double bad_ = 1.5, double cutoff_ = 8) :
  Term(cutoff_), good(good_), bad(bad_), type(TermType::ChargeIndependent) {
    term_name = "hydrophobic(g=" + to_string(good) + ",_b=" + to_string(bad)
            + ",_c=" + to_string(cutoff) + ")";
  }

  ~Hydrophobic() {
  }

  TermType get_term_type() {
    return type;
  }

  double eval(smt t1, smt t2, double r) {
    if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
      return slope_step(bad, good, r - optimal_distance(t1, t2));
    else
      return 0;
  }

  double eval_deriv(const smt t1, const smt t2, double r) {
    if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
      return slope_step_deriv(bad, good, r - optimal_distance(t1, t2));
    else
      return 0;
  }

  Term *copy() {
    return new Hydrophobic(*this);
  }

};

class Non_dir_h_bond : public Term {
public:
  TermType type;
  double good;
  double bad;

  Non_dir_h_bond(double good_ = -0.7, double bad_ = 0, double cutoff_ = 8) :
  Term(cutoff_), good(good_), bad(bad_), type(TermType::ChargeIndependent) {
    term_name = std::string("non_dir_h_bond(g=") + to_string(good) + ",_b="
            + to_string(bad) + ",_c=" + to_string(cutoff) + ")";
  }

  ~Non_dir_h_bond() {
  }

  TermType get_term_type() {
    return type;
  }

  double eval(smt t1, smt t2, double r) {
    if (xs_h_bond_possible(t1, t2))
      return slope_step(bad, good, r - optimal_distance(t1, t2));
    return 0;
  }

  double eval_deriv(const smt t1, const smt t2, double r) {
    if (xs_h_bond_possible(t1, t2))
      return slope_step_deriv(bad, good, r - optimal_distance(t1, t2));
    return 0;
  }

  Term *copy() {
    return new Non_dir_h_bond(*this);
  }
};

class Num_tors_div : public Term {
public:
  TermType type;

  Num_tors_div() : type(TermType::ConfIndependent) {
    term_name = "num_tors_div";
  }

  ~Num_tors_div() {
  }

  TermType get_term_type() {
    return type;
  }

  double eval(ConfIndependentInputs &in, double x, double w) {
    double new_w = 0.1 * (w + 1); // w is in [0..0.2]
    return smooth_div(x, 1 + new_w * in.num_tors / 5.0);
  }

  Term *copy() {
    return new Num_tors_div(*this);
  }
};
