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

#include "bond_order_determination.h"
#include "graph/cycle_enumerator.h"
#include "graph/graph.h"

using namespace std;
using namespace hess;

bool has_nonsp2(vector <int> & atoms, hess::Molecule* mol) {
  for (int i = 0; i < atoms.size(); i++) {
    if (mol->get_atom(atoms[i])->hybtype != 2)
      return true;
  }
  return false;
}

bool has_aromatic_bond(vector <int> & rings_bonds, hess::Molecule* mol) {
  for (int i = 0; i < rings_bonds.size(); i++) {
    if (mol->getBondOrder(rings_bonds[i]) == 4) {
      return true;
    }
  }
  return false;
}

void determine_bond_orders_in_aromatic_rings(hess::Molecule* mol, vector <vector <int>>&rings_bonds, vector <vector <int>>&rings_atoms) {
  for (int i = 0; i < rings_bonds.size(); i++) {
    for (int j = 0; j < rings_atoms[i].size(); j++) {
      int a_id = rings_atoms[i][j];
      int b_id = rings_atoms[i][(j + 1) % rings_atoms[i].size()];
      int bond_id = mol->findEdgeIndex(a_id, b_id);
      mol->get_bond(bond_id)->in_ring = true;
    }
  }
  int i = 0;
  while (i != rings_bonds.size()) {
    if (has_nonsp2(rings_atoms[i], mol)) {
      rings_bonds.erase(rings_bonds.begin() + i);
      rings_atoms.erase(rings_atoms.begin() + i);
    } else {
      i++;
    }
  }
  for (vector<int>& ring_bonds : rings_bonds) {
    for (int bond : ring_bonds) {
      mol->setBondOrder(bond, 4);
    }
  }

  indigo::AromaticityOptions arom_options(indigo::AromaticityOptions::GENERIC);

  mol->dearomatize(arom_options);
  for (int i = 0; i < rings_bonds.size(); i++) {
    if (has_aromatic_bond(rings_bonds[i], mol)) {
      for (int j = 0; j < rings_bonds[i].size(); j++) {
        if (mol->getBondOrder(rings_bonds[i][j]) == 4) {
          mol->setBondOrder(rings_bonds[i][j], 1);
        }
      }
      continue;
    }

    for (int j = 0; j < rings_atoms[i].size(); j++)
      mol->get_atom(rings_atoms[i][j])->aromatic = true;
  }

  for (int i = mol->vertexBegin(); i != mol->vertexEnd(); i = mol->vertexNext(i)) {
    int expl_val = 0;

    const indigo::Vertex &v = mol->getVertex(i);
    for (int j = v.neiBegin(); j != v.neiEnd(); j = v.neiNext(j)) {
      int order = mol->getBondOrder(v.neiEdge(j));

      if (order == indigo::BOND_AROMATIC)
        ; // should not happen
      else
        expl_val += order;
    }
    mol->atoms[i].valence = expl_val;
  }
}

double get_torsion_angle(int a_id, int b_id, int c_id, int d_id, hess::Molecule* mol) {
  hess::Atom *a = mol->get_atom(a_id);
  hess::Atom *b = mol->get_atom(b_id);
  hess::Atom *c = mol->get_atom(c_id);
  hess::Atom *d = mol->get_atom(d_id);
  hess::Vec3d a_v = hess::Vec3d(a->x, a->y, a->z);
  hess::Vec3d b_v = hess::Vec3d(b->x, b->y, b->z);
  hess::Vec3d c_v = hess::Vec3d(c->x, c->y, c->z);
  hess::Vec3d d_v = hess::Vec3d(d->x, d->y, d->z);
  double torsion;
  hess::Vec3d b1, b2, b3, c1, c2, c3, b2xb3, b1xb2;
  b1.diff(a_v, b_v);
  b2.diff(b_v, c_v);
  b3.diff(c_v, d_v);
  c1.cross(b1, b2);
  c2.cross(b2, b3);
  c3.cross(c1, c2);
  if (c1.length() * c2.length() < 0.001) {
    torsion = 0.0;
    return torsion;
  }
  double rb2 = sqrt(hess::Vec3d::dot(b2, b2));
  b2xb3.cross(b2, b3);
  b1xb2.cross(b1, b2);
  hess::Vec3d buf(0.0, 0.0, 0.0);
  buf.scaled(b1, rb2);
  torsion = -atan2(hess::Vec3d::dot(buf, b2xb3), hess::Vec3d::dot(b1xb2, b2xb3));
  return (torsion * RAD_TO_DEG);
}

void check_torsion_angles(hess::Molecule* mol, vector <vector <int>>&rings_atoms) {
  double torsions = 0.0;
  for (int i = 0; i < rings_atoms.size(); i++) {
    if (rings_atoms[i].size() == 5) {
      int a_0 = rings_atoms[i][0];
      int a_1 = rings_atoms[i][1];
      int a_2 = rings_atoms[i][2];
      int a_3 = rings_atoms[i][3];
      int a_4 = rings_atoms[i][4];
      vector <int> five_atoms = {a_0, a_1, a_2, a_3, a_4};
      double ts_1 = fabs(get_torsion_angle(five_atoms[0], five_atoms[1], five_atoms[2], five_atoms[3], mol));
      double ts_2 = fabs(get_torsion_angle(five_atoms[1], five_atoms[2], five_atoms[3], five_atoms[4], mol));
      double ts_3 = fabs(get_torsion_angle(five_atoms[2], five_atoms[3], five_atoms[4], five_atoms[0], mol));
      double ts_4 = fabs(get_torsion_angle(five_atoms[3], five_atoms[4], five_atoms[0], five_atoms[1], mol));
      double ts_5 = fabs(get_torsion_angle(five_atoms[4], five_atoms[0], five_atoms[1], five_atoms[2], mol));
      torsions = (ts_1 + ts_2 + ts_3 + ts_4 + ts_5) / 5.0;
      if (torsions <= 7.5) {
        for (int atom_id : five_atoms) {
          hess::Atom *atom = mol->get_atom(atom_id);
          if (mol->getVertex(atom_id).degree() == 2) {
            atom->hybtype = 2;
          }
        }
      }
    } else if (rings_atoms[i].size() == 6) {
      int a_0 = rings_atoms[i][0];
      int a_1 = rings_atoms[i][1];
      int a_2 = rings_atoms[i][2];
      int a_3 = rings_atoms[i][3];
      int a_4 = rings_atoms[i][4];
      int a_5 = rings_atoms[i][5];
      vector <int> six_atoms = {a_0, a_1, a_2, a_3, a_4, a_5};
      double ts_1 = fabs(get_torsion_angle(six_atoms[0], six_atoms[1], six_atoms[2], six_atoms[3], mol));
      double ts_2 = fabs(get_torsion_angle(six_atoms[1], six_atoms[2], six_atoms[3], six_atoms[4], mol));
      double ts_3 = fabs(get_torsion_angle(six_atoms[2], six_atoms[3], six_atoms[4], six_atoms[5], mol));
      double ts_4 = fabs(get_torsion_angle(six_atoms[3], six_atoms[4], six_atoms[5], six_atoms[0], mol));
      double ts_5 = fabs(get_torsion_angle(six_atoms[4], six_atoms[5], six_atoms[0], six_atoms[1], mol));
      double ts_6 = fabs(get_torsion_angle(six_atoms[5], six_atoms[0], six_atoms[1], six_atoms[2], mol));
      torsions = (ts_1 + ts_2 + ts_3 + ts_4 + ts_5 + ts_6) / 6.0;
      if (torsions <= 12) {
        for (int atom_id : six_atoms) {
          hess::Atom *atom = mol->get_atom(atom_id);
          if (mol->getVertex(atom_id).degree() == 2 || mol->getVertex(atom_id).degree() == 3) {
            atom->hybtype = 2;
          }
        }
      }
    }
  }
}

struct _Data {

  _Data(vector <vector <int>> &a, vector <vector <int>> &b) : rings_atoms(a), rings_bonds(b) {
  }
  vector <vector <int>> &rings_atoms;
  vector <vector <int>> &rings_bonds;

};

bool _handleCycle(indigo::Graph& graph, const indigo::Array<int>& vertices, const indigo::Array<int>& edges, void* context) {
  _Data* data = (_Data*) context;

  std::vector<int> v;
  for (int i = 0; i < vertices.size(); i++)
    v.push_back(vertices[i]);
  data->rings_atoms.push_back(move(v));

  std::vector<int> e;
  for (int i = 0; i < edges.size(); i++)
    e.push_back(edges[i]);
  data->rings_bonds.push_back(move(e));

  return true;
}

void extract_rings(vector <vector <int>> &rings_bonds, vector <vector <int>> &rings_atoms, hess::Molecule *mol) {
  indigo::CycleEnumerator enumerator(*mol);
  _Data data(rings_atoms, rings_bonds);
  enumerator.min_length = 1;
  enumerator.max_length = 100;
  enumerator.context = &data;
  enumerator.cb_handle_cycle = _handleCycle;
  enumerator.process();
}
