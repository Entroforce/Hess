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

#include "global_optimizer.h"
#include "transformation.h"
#include "logger.h"

using namespace std;

void moveProtein(hess::Molecule* ord_rec, double* box) {
  double xc = box[0], yc = box[1], zc = box[2];
  for (int a_id = ord_rec->vertexBegin(); a_id != ord_rec->vertexEnd(); a_id = ord_rec->vertexNext(a_id)) {
    hess::Atom* a = ord_rec->get_atom(a_id);
    a->x -= xc;
    a->y -= yc;
    a->z -= zc;
  }
}

void calcLigandRootCenter(double* ligand_center, hess::Molecule * ligand) {
  simplified_tree tr(ligand->atoms.size());
  tree_init(tr, ligand);
  int count = 0;
  hess::Vec3d average;
  calc_center_mass_frag(ligand, tr, tr.rootid, count, average);
  average.x /= count;
  average.y /= count;
  average.z /= count;
  ligand_center[0] = average.x;
  ligand_center[1] = average.y;
  ligand_center[2] = average.z;
}

void moveLigandToCenter(double* ligand_center, hess::Molecule* ligand) {
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    hess::Atom* a = ligand->get_atom(a_id);
    a->x -= ligand_center[0];
    a->y -= ligand_center[1];
    a->z -= ligand_center[2];
    a->x_orign -= ligand_center[0];
    a->y_orign -= ligand_center[1];
    a->z_orign -= ligand_center[2];
  }
}

void delete_nonpolar_hydrogens(hess::Molecule *mol) {
  vector<int> del_id;
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    if (mol->atom_is_nonpolar_hydrogen(a_id) && is_suppressible_hydrogen(mol, a_id)) {
      del_id.push_back(a_id);
    }
  }
  for (int a_id : del_id) {
    mol->remove_atom(a_id);
  }
}

void assign_bond_rots(hess::Molecule *mol) {
  for (int bond_id = mol->edgeBegin(); bond_id != mol->edgeEnd(); bond_id = mol->edgeNext(bond_id)) {
    hess::Bond* bond = mol->get_bond(bond_id);
    if (is_rot_bond(mol, bond_id)) {
      bond->rotatable = true;
    }
  }
}

Optimizable_molecule::Optimizable_molecule(hess::Molecule* _lig, hess::Molecule* _prot, double* box, const char* _optimize, double _granularity, unsigned _seed) :
grid_new(num_atom_types()), grid_deriv_new(num_atom_types()) {
  scoring = "vinardo";

  assign_atom_types(_lig, string(scoring));
  assign_atom_types(_prot, string(scoring));
  adjust_atom_types(_lig);
  adjust_atom_types(_prot);
  delete_nonpolar_hydrogens(_lig);
  delete_nonpolar_hydrogens(_prot);
  assign_bond_rots(_lig);

  ligand = new hess::Molecule(*((hess::Molecule*)_lig));
  receptor = new hess::Molecule(*((hess::Molecule*)_prot));
  moveProtein(receptor, box);
  double ligand_center[3] = {0};
  calcLigandRootCenter(ligand_center, ligand);
  moveLigandToCenter(ligand_center, ligand);
  find_ligand_pairs(ligand);
  set_insubtree(ligand);  
  set_parents_to_tree(ligand);

  in = ConfIndependentInputs(ligand);
  type = 0;
  optimize = _optimize;
  ligand->coords.resize(ligand->get_atoms_count());
  rot_bonds_count = ligand->get_rotatable_bonds_count();
  solve_count = rot_bonds_count + 6;
  fixedcon.resize(rot_bonds_count + 6);
  seed = _seed;
  simplified_tree _tr(ligand->atoms.size());
  tr = _tr;
  tree_init(tr, ligand);
  for (int i = 0; i < tr.bond_to_parent.size(); i++) {
    if (tr.bond_to_parent[i]) {
      rotsmap.push_back(i);
    }
  }
  for (int i = 0; i < tr.sz; i++) {
    if (tr.bond_to_parent[i])
      encoding.push_back(i);
  }
  for (int i = tr.sz; i < tr.sz + 6; i++) {
    encoding.push_back(i);
  }
  encoding_inv.resize(ligand->get_atoms_count() + 6);
  for (int i = 0; i < encoding.size(); i++) {
    encoding_inv[encoding[i]] = i;
  }
  size_x = box[3];
  size_y = box[4];
  size_z = box[5];
  xc = box[0];
  yc = box[1];
  zc = box[2];
  this->ligand_center = ligand_center;
  granularity = _granularity;
  const double span[3] = {this->size_x, this->size_y, this->size_z};
  for (int i = 0; i < 3; i++) {
    this->gd[i][2] = size_t(std::ceil(span[i] / granularity)); //n
    double real_span = granularity * this->gd[i][2];
    this->gd[i][0] = this->ligand_center[i] - real_span / 2; // begin
    this->corner1[i] = this->gd[i][0];
    this->gd[i][1] = this->gd[i][0] + real_span; // end
    this->corner2[i] = this->gd[i][1];
  }
  fprintf(hess::Logger::get_instance().get_stream(), "Grid granularity: %.4g\n", granularity);
  dim_x = gd[0][2] + 1;
  dim_y = gd[1][2] + 1;
  dim_z = gd[2][2] + 1;
  dim_fl_minus_1[0] = dim_x - 1;
  dim_fl_minus_1[1] = dim_y - 1;
  dim_fl_minus_1[2] = dim_z - 1;
  for (int i = 0; i < 3; i++) {
    init[i] = gd[i][0];
    range[i] = gd[i][1] - gd[i][0];
    factor[i] = dim_fl_minus_1[i] / range[i];
    factor_inv[i] = 1 / factor[i];
  }
  init_vec.set(init[0], init[1], init[2]);
  factor_vec.set(factor[0], factor[1], factor[2]);
  factor_inv_vec.set(factor_inv[0], factor_inv[1], factor_inv[2]);
  for (int a_id = this->ligand->vertexBegin(); a_id != this->ligand->vertexEnd(); a_id = this->ligand->vertexNext(a_id)) {
    int n = num_atom_types();
    hess::Atom* a = ligand->get_atom(a_id);
    smt t = a->atom_type;
    if (t < n && ligand_types.find(t) == ligand_types.end() && !is_hydrogen(t)) {
      ligand_types.insert(t);
      ligand_types_vec.push_back(t);
      grid_new[t].resize(dim_x, dim_y, dim_z);
      grid_deriv_new[t].resize(dim_x, dim_y, dim_z);
    }
  }
  int heuristic = ligand->get_atoms_count() + 10 * (rot_bonds_count + 6);
  num_steps = unsigned(70 * 3 * (50 + heuristic) / 2);
}

Optimizable_molecule::~Optimizable_molecule() {
  delete this->ligand;
  delete this->receptor;
  encoding.clear();
  encoding_inv.clear();
  rotsmap.clear();
  grid_new.clear();
  grid_deriv_new.clear();
  ligand_types_vec.clear();
  for (int i = 0; i < result_mols.size(); i++) {
    delete result_mols[i];
  }
  result_mols.clear();

}

double Optimizable_molecule::calc_derivative_fast_type_0_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
  //    fast_grad_comp_count += 1;
  double res_energy = fast_gradient(vinardo_sc, this, x, grad, encoding_inv);
  return res_energy;
}

double Optimizable_molecule::calc_derivative_fast_type_1_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
  //    fast_grad_comp_count += 1;
  Eigen::VectorXd x0 = Eigen::VectorXd::Zero(fixedcon.size());
  for (int i = 0; i < fixedcon.size(); i++) {
    x0[i] = fixedcon[i];
  }
  for (int i = 0; i < 6; i++) {
    x0[this->rot_bonds_count + i] = x[i];
  }
  double res_energy = fast_gradient(vinardo_sc, this, x0, grad, encoding_inv);
  return res_energy;
}

double Optimizable_molecule::calc_derivative_fast_upd(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
  double energy = 0;
  if (!this->type)
    energy = this->calc_derivative_fast_type_0_upd(x, grad);
  else if (this->type == 1)
    energy = this->calc_derivative_fast_type_1_upd(x, grad);
  return energy;
}

double Optimizable_molecule::operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad) {
  for (int i = 0; i < x.size(); i++) {
    if (isnan(x[i])) {
      for (int i = 0; i < grad.size(); i++) {
        grad[i] = slope;
      }
      return slope;
    }
  }
  for (int i = 0; i < grad.size(); i++) {
    grad[i] = 0;
  }
  double energy = this->calc_derivative_fast_upd(x, grad);
  return energy;
}

void Optimizable_molecule::set_max_bfgs_iterations(int max_bfgs_iterations) {
  this->max_bfgs_iterations = max_bfgs_iterations;
}

void Optimizable_molecule::fill_grid() {
  printf("Building grid...\n");
  const int n = num_atom_types();
  double max_cutoff_sqr = vinardo_sc.get_max_cutoff_sqr();
  std::vector<hess::Atom*> prot_atoms;
  for (int b_id = receptor->vertexBegin(); b_id != receptor->vertexEnd(); b_id = receptor->vertexNext(b_id))
    prot_atoms.push_back(receptor->get_atom(b_id));
  double ax, ay, az, r2;
  double all_t_acc[specific_atom_type::type::NumTypes] = {0};
  for (int x = 0; x < this->dim_x; x++) {
    ax = init[0] + factor_inv[0] * x;
    for (int y = 0; y < this->dim_y; y++) {
      ay = init[1] + factor_inv[1] * y;
      for (int z = 0; z < this->dim_z; z++) {
        az = init[2] + factor_inv[2] * z;
        for (hess::Atom* b : prot_atoms) {
          smt t2 = b->atom_type;
          if (t2 >= n || is_hydrogen(t2)) {
            continue;
          }
          r2 = (b->x - ax) * (b->x - ax) + (b->y - ay) * (b->y - ay) + (b->z - az) * (b->z - az);
          if (r2 < max_cutoff_sqr)
            for (smt t : this->ligand_types_vec)
              all_t_acc[t] += vinardo_sc.eval_fast(t, t2, sqrt(r2));
        }
        for (smt t : this->ligand_types_vec) {
          this->grid_new[t](x, y, z) = all_t_acc[t];
          all_t_acc[t] = 0;
        }

      }
    }
  }
}

void Optimizable_molecule::fill_grid_with_derivs() {
  fprintf(hess::Logger::get_instance().get_stream(), "Building grid with gradients...\n");
  const int n = num_atom_types();
  double max_cutoff_sqr = vinardo_sc.get_max_cutoff_sqr();
  std::vector<hess::Atom*> prot_atoms;
  for (int b_id = receptor->vertexBegin(); b_id != receptor->vertexEnd(); b_id = receptor->vertexNext(b_id))
    prot_atoms.push_back(receptor->get_atom(b_id));
  double ax, ay, az, r2;
  double all_t_acc[specific_atom_type::type::NumTypes] = {0};
  hess::Vec3d acc_d_terms[specific_atom_type::type::NumTypes] = {
    {0.0, 0.0, 0.0}};
  for (int x = 0; x < this->dim_x; x++) {
    ax = init[0] + factor_inv[0] * x;
    for (int y = 0; y < this->dim_y; y++) {
      ay = init[1] + factor_inv[1] * y;
      for (int z = 0; z < this->dim_z; z++) {
        az = init[2] + factor_inv[2] * z;
        for (hess::Atom* b : prot_atoms) {
          hess::Vec3d b_vec = b->get_vector();
          smt t2 = b->atom_type;
          if (t2 >= n || is_hydrogen(t2)) {
            continue;
          }
          hess::Vec3d approx_vector(ax, ay, az);
          r2 = (b->x - ax) * (b->x - ax) + (b->y - ay) * (b->y - ay) + (b->z - az) * (b->z - az);
          if (r2 < max_cutoff_sqr)
            for (smt t : this->ligand_types_vec) {
              all_t_acc[t] += vinardo_sc.eval_fast(t, t2, sqrt(r2));
              hess::Vec3d d_dterms_mul = vinardo_sc.calc_term_deriv_part_dist(approx_vector, b_vec);
              double d_terms = vinardo_sc.calc_term_deriv_part_terms(t, t2, sqrt(r2));
              d_dterms_mul.x *= d_terms;
              d_dterms_mul.y *= d_terms;
              d_dterms_mul.z *= d_terms;
              acc_d_terms[t].add(d_dterms_mul);

            }
        }
        for (smt t : this->ligand_types_vec) {
          this->grid_new[t](x, y, z) = all_t_acc[t];
          this->grid_deriv_new[t](x, y, z) = acc_d_terms[t];
          all_t_acc[t] = 0;
          acc_d_terms[t] = {0, 0, 0};
        }

      }
    }
  }
}

void Optimizable_molecule::assign_hunt_cap() {
  this->cap[0] = this->hunt_cap[0];
  this->cap[1] = this->hunt_cap[1];
}

void Optimizable_molecule::assign_aut_cap() {
  this->cap[0] = this->aut_cap[0];
  this->cap[1] = this->aut_cap[1];
}
