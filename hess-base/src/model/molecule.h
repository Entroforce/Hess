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

#include "exception.h"
#include "vector3.h"
#include "molecule/molecule.h"
#include "scoring_function/constants.h"

using std::vector;
using namespace std;

typedef const indigo::Vertex& ver;
namespace hess {

  class Atom {
  public:
    Atom(double x = 0.0, double y = 0.0, double z = 0.0);
    string get_type();
    hess::Vec3d get_vector();

    int valence; // explicit valence
    double x;
    double y;
    double z;
    double x_orign; //_orign for optimizations
    double y_orign;
    double z_orign;
    int hybtype;
    bool in_ring;
    bool aromatic;
    smt atom_type;
    int num;
  };

  class Bond {
  public:
    Bond(double length = 0.0);
    double length;
    bool rotatable;
    bool in_ring;
  };

  class Molecule : public indigo::Molecule {
  public:
    Molecule();
    ~Molecule();
    Molecule(Molecule &mol);
    int get_atoms_count() const;
    void calc_coords();
    const vector<pair<int, int>>&get_lig_pairs() const;
    const vector<vector<int>>&get_insubt() const;
    void set_insubt(vector<vector<int>>&& ins);
    void set_tree_parents(vector<int>&& tree_parents);
    const vector<int> &get_tree_parents() const;
    void charge_calculate();
    int bonds_count() const;
    Atom* get_atom(int id);
    const Atom* get_atom(int id)const;
    Bond* get_bond(int id);
    const Bond* get_bond(int id) const;
    int add_atom(double x = 0.0, double y = 0.0, double z = 0.0, const int element_num = 0);
    int add_bond(int beg, int end, int order, double length);
    int bond_atom_beg(int bond_id) const;
    int bond_atom_end(int bond_id) const;
    void remove_atom(int id);
    void remove_bond(int id);
    bool atom_is_nonpolar_hydrogen(int atom_id);
    void set_rings(const vector<vector<int>> &rings);
    vector <vector <int>> get_rings();
    void set_aromatic_rings(vector<vector<int>>&& aromatic_rings);
    vector <vector <int>> get_aromatic_rings();
    int get_rotatable_bonds_count();


    int determine_order(int ia, int ib, double* length_bound = nullptr);
    double get_correct_bond_rad(int idx);

    double pH;
    int total_charge;
    vector<Atom> atoms;
    vector<Bond> bonds;
    vector <vector <int>> rings;
    vector <vector <int>> aromatic_rings;
    vector<pair<int, int>> pairs;
    vector<vector<int>> insubt;
    vector<int> tree_parents;
    vector<hess::Vec3d> coords;
    char **txtatm;
  };
} // namespace hess
