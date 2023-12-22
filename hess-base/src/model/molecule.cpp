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

#include "model/molecule.h"
#include "scoring_function/features.h"
#include "hess.h"
#include "elements.h"

using namespace hess;

hess::Atom::Atom(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->x_orign = x;
  this->y_orign = y;
  this->z_orign = z;
  this->valence = 0;
  this->in_ring = false;
  this->hybtype = 0;
  this->aromatic = false;
}

hess::Vec3d hess::Atom::get_vector() {
  return hess::Vec3d(this->x, this->y, this->z);
}

namespace specific_atom_type {
  info data[NumTypes] = {
    {},};
}

Bond::Bond(double length) {
  this->length = length;
  this->rotatable = false;
  this->in_ring = false;
}

Molecule::Molecule() {
  this->pH = -1.0;
  total_charge = 0;
};

Molecule::~Molecule() {
}

void Molecule::charge_calculate() {
  int total = 0;
  for (int a_id = this->vertexBegin(); a_id != this->vertexEnd(); a_id = this->vertexNext(a_id)) {
    total += getAtomCharge(a_id);
  }
  this->total_charge = total;
}

Molecule::Molecule(Molecule &mol) {
  this->pH = mol.pH;
  this->total_charge = mol.total_charge;

  map<int, int> map_id;
  for (int a_id = mol.vertexBegin(); a_id != mol.vertexEnd(); a_id = mol.vertexNext(a_id)) {
    int copy_id = this->add_atom(mol.get_atom(a_id)->x, mol.get_atom(a_id)->y, mol.get_atom(a_id)->z, mol.getAtomNumber(a_id));
    atoms[copy_id] = mol.atoms[a_id];
    map_id[a_id] = copy_id;
  }

  for (int bond_id = mol.edgeBegin(); bond_id != mol.edgeEnd(); bond_id = mol.edgeNext(bond_id)) {
    int a_id = mol.bond_atom_beg(bond_id);
    int b_id = mol.bond_atom_end(bond_id);
    int curr_a_id = map_id[a_id];
    int curr_b_id = map_id[b_id];
    int copy_id = this->add_bond(curr_a_id, curr_b_id, mol.getBondOrder(bond_id), mol.bonds[bond_id].length);
    bonds[copy_id] = mol.bonds[bond_id];
  }
  this->pairs = mol.pairs;
  this->insubt = mol.insubt;
  this->tree_parents = mol.tree_parents;
  this->coords = mol.coords;
  // ring orders are not recalculated!!!
  for (vector<int>& ring : mol.rings) {
    vector<int> new_ring;
    for (int old_id : ring) {
      new_ring.push_back(map_id[old_id]);
    }
    this->rings.push_back(new_ring);
  }
  //also that aromatic rings
  for (vector<int>& aromatic_ring : mol.aromatic_rings) {
    vector<int> new_aromatic_ring;
    for (int old_id : aromatic_ring) {
      new_aromatic_ring.push_back(map_id[old_id]);
    }
    this->aromatic_rings.push_back(new_aromatic_ring);
  }
}

bool Molecule::atom_is_nonpolar_hydrogen(int atom_id) {
  if (getAtomNumber(atom_id) != 1)
    return false;

  ver v = this->getVertex(atom_id);

  for (int id = v.neiBegin(); id != v.neiEnd(); id = v.neiNext(id))
    if (getAtomNumber(v.neiVertex(id)) == 6)
      return true;

  return false;
}

int Molecule::get_atoms_count() const {
  return this->vertexCount();
}

void Molecule::calc_coords() {
  for (int a_id = this->vertexBegin(); a_id != this->vertexEnd(); a_id = this->vertexNext(a_id)) {
    Atom *a = this->get_atom(a_id);
    coords[a_id] = a->get_vector();
  }
}

const vector<pair<int, int>>&Molecule::get_lig_pairs() const {
  return this->pairs;
}

const vector<vector<int>>&Molecule::get_insubt() const {
  return this->insubt;
}

void Molecule::set_insubt(vector<vector<int>>&& ins) {
  this->insubt = ins;
}

void Molecule::set_tree_parents(vector<int>&& tree_parents) {
  this->tree_parents = tree_parents;
}

const vector<int>& Molecule::get_tree_parents() const {
  return this->tree_parents;
}

int Molecule::bonds_count() const {
  return this->edgeCount();
}

Atom* Molecule::get_atom(int id) {
  if (id >= this->atoms.size()) {
    throw HessException("Atom index out of bounds");
  } else {
    return &this->atoms[id];
  }
}

const Atom* Molecule::get_atom(int id) const {
  return get_atom(id);
}

Bond* Molecule::get_bond(int id) {
  if (id >= this->bonds.size()) {
    throw HessException("Bond index out of bounds");
  } else {
    return &this->bonds[id];
  }
}

const Bond* Molecule::get_bond(int id) const {
  return get_bond(id);
}

int Molecule::add_atom(double x, double y, double z, const int element_num) {
  int id = addAtom(element_num);

  if (id >= atoms.size())
    atoms.resize(id + 1);

  memset(&atoms[id], 0, sizeof (Atom));

  atoms[id].x = x;
  atoms[id].y = y;
  atoms[id].z = z;
  atoms[id].x_orign = x;
  atoms[id].y_orign = y;
  atoms[id].z_orign = z;

  setAtomXyz(id, x, y, z);
  return id;
}

int Molecule::add_bond(int beg, int end, int order, double length) {
  int id = addBond(beg, end, order);

  if (id >= bonds.size())
    bonds.resize(id + 1);

  memset(&bonds[id], 0, sizeof (Bond));
  bonds[id].length = length;
  return id;
}

int Molecule::bond_atom_beg(int bond_id) const {
  return this->getEdge(bond_id).beg;
}

int Molecule::bond_atom_end(int bond_id) const {
  return this->getEdge(bond_id).end;
}

void Molecule::remove_atom(int id) {
  ver atom_ver = getVertex(id);
  for (int idd = atom_ver.neiBegin(); idd != atom_ver.neiEnd(); idd = atom_ver.neiNext(idd)) {
    int nei_id = atom_ver.neiVertex(idd);
    int bond_id = this->findEdgeIndex(id, nei_id);
    Atom *nei = this->get_atom(nei_id);
    nei->valence--;
    if (getBondOrder(bond_id) == 2) {
      nei->valence--;
    } else if (getBondOrder(bond_id) == 3) {
      nei->valence--;
      nei->valence--;
    }
  }
  removeAtom(id);
}

void Molecule::remove_bond(int id) {
  int a_id = this->bond_atom_beg(id);
  int b_id = this->bond_atom_end(id);
  Atom* a = get_atom(a_id);
  Atom* b = get_atom(b_id);
  a->valence--;
  b->valence--;
  if (getBondOrder(id) == 2) {
    a->valence--;
    b->valence--;
  } else if (getBondOrder(id) == 3) {
    a->valence--;
    b->valence--;
    a->valence--;
    b->valence--;
  }

  removeBond(id);
}

bool Molecule::atom_is_nonpolar_hydrogen(int atom_id);

void Molecule::set_rings(const vector<vector<int>> &rings) {
  this->rings = rings;
}

vector <vector <int>> Molecule::get_rings() {
  return this->rings;
}

void Molecule::set_aromatic_rings(vector<vector<int>>&& aromatic_rings) {
  this->aromatic_rings = aromatic_rings;
}

vector <vector <int>> Molecule::get_aromatic_rings() {
  return this->aromatic_rings;
}

int Molecule::get_rotatable_bonds_count() {
  int count = 0;
  for (int bond_id = this->edgeBegin(); bond_id != this->edgeEnd(); bond_id = this->edgeNext(bond_id)) {
    const Bond* bond = this->get_bond(bond_id);
    if (bond->rotatable)
      count += 1;
  }
  return count;
}

int Molecule::determine_order(int ia, int ib, double* length_bound) {
  double length = distance(atoms[ia], atoms[ib]);
  if (length_bound != nullptr)
    *length_bound = length;
  double up_lim = 0.45;
  double down_lim = 0.4;
  int order = 0;
  double cov_a = GetCovalentRad(getAtomNumber(ia));
  double cov_b = GetCovalentRad(getAtomNumber(ib));

  if (cov_a == 0.0 || cov_b == 0.0) {
    fprintf(hessGetStream(), "Did not find rad in obabel table!\n");
  }
  if (length < cov_a + cov_b + up_lim && length > cov_a + cov_b - down_lim)
    order = 1;

  return order;
}

double Molecule::get_correct_bond_rad(int idx) {
  int atom_num = getAtomNumber(idx);
  int type = atoms[idx].hybtype;
  if (atom_num < 0)
    return 1.0; // TODO: ???? can not happen
  double rad = GetCovalentRad(atom_num);
  if (type == 2)
    return rad * 0.95;
  else if (type == 1)
    return rad * 0.90;
  return rad;
}
