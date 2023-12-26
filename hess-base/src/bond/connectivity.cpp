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
#include "molecule/elements.h"
#include "elements.h"
#include "hess.h"

using namespace std;
using namespace indigo;
using namespace hess;

int count_explicit_h(int a_id, hess::Molecule* mol) {
  int h_count = 0;
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id))
    if (mol->getAtomNumber(atom_ver.neiVertex(id)) == 1)
      h_count++;
  return h_count;
}

int count_free_oxygens(hess::Molecule* mol, const int a_id) {
  int count = 0;
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (mol->getAtomNumber(nei_id) == ELEM_O && nonhydrobond_count(mol, nei_id) == 1)
      count++;

  }
  return count;
}


//! \return The number of sulfur atoms connected that only have one heavy valence

int get_count_free_sulfurs(hess::Molecule* mol, int a_id) {
  int count = 0;
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (mol->getAtomNumber(nei_id) == ELEM_S && nonhydrobond_count(mol, nei_id) == 1)
      count++;
  }
  return count;
}

int nonhydrobond_count(hess::Molecule* mol, int a_id) {
  return mol->getVertex(a_id).degree() - count_explicit_h(a_id, mol);
}

bool is_carboxyl_oxygen(hess::Molecule* mol, int a_id) {
  if (mol->getAtomNumber(a_id) != ELEM_O) {
    return false;
  }
  if (nonhydrobond_count(mol, a_id) != 1)
    return false;


  int carb_id = -1;
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (mol->getAtomNumber(nei_id) == ELEM_C) {
      carb_id = nei_id;
      break;
    }
  }
  if (carb_id == -1)
    return false;

  int count_free_oxygs = count_free_oxygens(mol, carb_id);
  int count_free_sulfs = get_count_free_sulfurs(mol, carb_id);
  if (!(count_free_oxygs == 2) && !(count_free_oxygs == 1 && count_free_sulfs == 1)) {
    return false;
  }
  return true;

}

//! \return returns the location of the added hydrogen

hess::Vec3d get_new_bonds_vector(hess::Molecule* mol, int a_id, double length) {
  Atom *a = mol->get_atom(a_id);
  const hess::Vec3d VX(1.0, 0.0, 0.0);
  hess::Vec3d v1, v2, newbond;
  hess::Vec3d bond1(0.0, 0.0, 0.0);
  hess::Vec3d bond2(0.0, 0.0, 0.0);
  hess::Vec3d bond3(0.0, 0.0, 0.0);
  hess::Vec3d bond4(0.0, 0.0, 0.0);
  hess::Vec3d bond5(0.0, 0.0, 0.0);
  if (mol->getVertex(a_id).degree() == 0) {
    hess::Vec3d buf(0.0, 0.0, 0.0);
    buf.scaled(VX, length);
    newbond.sum(a->get_vector(), buf);
    return newbond;
  }
  if (mol->getVertex(a_id).degree() == 1) {
    bool is_carboxylate_O = is_carboxyl_oxygen(mol, a_id);
    ver atom_ver = mol->getVertex(a_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      Atom *nbr = mol->get_atom(nei_id);
      hess::Vec3d a_coord = a->get_vector();
      hess::Vec3d nbr_coord = nbr->get_vector();
      bond1.diff(a_coord, nbr_coord);
      if (nbr->hybtype == 1) {
        continue;
      }
      ver atom_ver2 = mol->getVertex(nei_id);
      for (int id2 = atom_ver2.neiBegin(); id2 != atom_ver2.neiEnd(); id2 = atom_ver2.neiNext(id2)) {
        int nei_id2 = atom_ver2.neiVertex(id2);
        Atom *nbr2 = mol->get_atom(nei_id2);
        if (nei_id2 != a_id) {
          bond2.diff(nbr->get_vector(), nbr2->get_vector());
          if (is_carboxylate_O && mol->getAtomNumber(nei_id2) == ELEM_O) {
            break;
          }
        }
      }
    }
    bond1.normalize();
    v1.cross(bond1, bond2);
    if (is_zero(bond2) || is_zero(v1)) {
      hess::Vec3d vrand(0.0, 0.0, 0.0);
      random_vector(vrand);
      double angle = fabs(acos(hess::Vec3d::dot(bond1, vrand)) * RAD_TO_DEG);
      while (angle < 45.0 || angle > 135.0) {
        random_vector(vrand);
        angle = fabs(acos(hess::Vec3d::dot(bond1, vrand)) * RAD_TO_DEG);
      }
      v1.cross(bond1, vrand);
    }
    v2.cross(bond1, v1);
    v2.normalize();
    if (a->hybtype == 1)
      newbond.set(bond1.x, bond1.y, bond1.z);
    else if (a->hybtype == 2) {
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.scaled(v2, tan(DEG_TO_RAD * 60.0));
      newbond.diff(bond1, buf);
    } else if (a->hybtype == 3) {
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.scaled(v2, tan(DEG_TO_RAD * (180.0 - 109.471)));
      newbond.diff(bond1, buf);
    } else if (a->hybtype == 4)
      newbond.set(bond1.x, bond1.y, bond1.z);
    else if (a->hybtype == 5) {
      newbond.set(bond1.x, bond1.y, bond1.z);
    }
    else if (a->hybtype == 6) {
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.scaled(v2, tan(DEG_TO_RAD * 90.0));
      newbond.diff(bond1, buf);
    }
    newbond.normalize();
    newbond.scale(length);
    newbond.add(a->get_vector());
    return newbond;
  }
  else if (mol->getVertex(a_id).degree() == 2) {
    ver atom_ver = mol->getVertex(a_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      Atom *nbr = mol->get_atom(nei_id);
      if (is_zero(bond1)) {
        bond1.diff(a->get_vector(), nbr->get_vector());
      } else if (is_zero(bond2)) {
        bond2.diff(a->get_vector(), nbr->get_vector());
      }
    }
    bond1.normalize();
    bond2.normalize();
    v1.sum(bond1, bond2);
    v1.normalize();
    if (a->hybtype == 2) {
      newbond.set(v1.x, v1.y, v1.z);
    } else if (a->hybtype == 3) {
      v2.cross(bond1, bond2);
      v2.normalize();
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.scaled(v2, tan(DEG_TO_RAD * (180.0 - 109.471)));
      newbond.diff(bond1, buf);
      hess::Vec3d buf1(0.0, 0.0, 0.0);
      buf1.scaled(v1, (sqrt(2.0) / 2.0));
      newbond.sum(v2, buf1);
    } else if (a->hybtype == 4 || a->hybtype == 5) {
      hess::Vec3d vrand(0.0, 0.0, 0.0);
      random_vector(vrand);
      double angle = fabs(acos(hess::Vec3d::dot(bond1, vrand)) * RAD_TO_DEG);
      while (angle < 45.0 || angle > 135.0) {
        random_vector(vrand);
        angle = fabs(acos(hess::Vec3d::dot(bond1, vrand)) * RAD_TO_DEG);
      }
      v1.cross(bond1, vrand);
      v1.normalize();
      newbond.set(v1.x, v1.y, v1.z);
    } else if (a->hybtype == 6) {
      v2.cross(bond1, bond2);
      newbond.set(v2.x, v2.y, v2.z);
    }
    newbond.normalize();
    newbond.scale(length);
    newbond.add(a->get_vector());
    return newbond;
  } else if (mol->getVertex(a_id).degree() == 3) {
    if (a->hybtype == 3) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        }
      }
      bond1.normalize();
      bond2.normalize();
      bond3.normalize();
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.sum(bond1, bond2);
      newbond.sum(bond3, buf);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
    if (a->hybtype == 4) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        }
      }
      bond3.normalize();
      newbond.set(bond3.x, bond3.y, bond3.z);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
    if (a->hybtype == 5) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        }
      }
      bond1.normalize();
      bond2.normalize();
      bond3.normalize();
      v1.cross(bond1, bond3);
      v1.normalize();
      hess::Vec3d buf(0.0, 0.0, 0.0);
      buf.scaled(bond3, tan(DEG_TO_RAD * 30.0));
      newbond.sum(v1, buf);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
    if (a->hybtype == 6) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        }
      }
      newbond.set(bond1.x, bond1.y, bond1.z);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
  } else if (mol->getVertex(a_id).degree() == 4) {
    if (a->hybtype == 6) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond4)) {
          bond4.diff(a->get_vector(), nbr->get_vector());
        }
      }
      newbond.set(bond2.x, bond2.y, bond2.z);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    } else if (a->hybtype == 5) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond4)) {
          bond4.diff(a->get_vector(), nbr->get_vector());
        }
      }
      bond1.normalize();
      bond2.normalize();
      bond3.normalize();
      bond4.normalize();
      v1.cross(bond1, bond3);
      v1.normalize();
      hess::Vec3d buf1(0.0, 0.0, 0.0);
      hess::Vec3d buf2(0.0, 0.0, 0.0);
      buf1.negation(v1);
      buf2.scaled(bond3, tan(DEG_TO_RAD * 30.0));
      newbond.sum(buf1, buf2);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
  }
  else if (mol->getVertex(a_id).degree() == 5) {
    if (a->hybtype == 6) {
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        Atom *nbr = mol->get_atom(nei_id);
        if (is_zero(bond1)) {
          bond1.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond2)) {
          bond2.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond3)) {
          bond3.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond4)) {
          bond4.diff(a->get_vector(), nbr->get_vector());
        } else if (is_zero(bond5)) {
          bond5.diff(a->get_vector(), nbr->get_vector());
        }
      }
      newbond.set(bond3.x, bond3.y, bond3.z);
      newbond.normalize();
      newbond.scale(length);
      newbond.add(a->get_vector());
      return newbond;
    }
  }
  random_vector(newbond);
  newbond.scale(length);
  newbond.add(a->get_vector());
  return newbond;
}

void unfold_hydrogens(hess::Molecule* mol) {
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    int h_count = mol->getImplicitH_NoThrow(a_id, 0);
    if (h_count != 0) {
      while (h_count > 0) {
        double length = mol->get_correct_bond_rad(a_id) + GetCovalentRad(1);
        hess::Vec3d h_coord = get_new_bonds_vector(mol, a_id, length);
        int h_id = mol->add_atom(h_coord.x, h_coord.y, h_coord.z, 1);
        hess::Atom *new_h = mol->get_atom(h_id);
        hess::Atom *a = mol->get_atom(a_id);
        a->valence++;
        new_h->valence++;
        mol->add_bond(a_id, h_id, 1, mol->get_correct_bond_rad(a_id) + GetCovalentRad(1));
        h_count--;
      }
    }
  }
}

//! \return true if the values ​​differ with an accuracy not less than the specified value.

bool is_aprox(double a, double b, double accuracy) {
  return abs(a - b) <= accuracy;
}

bool sort_atom_z(const pair<int, double> &a, const pair<int, double> &b) {
  return (a.second < b.second);
}


//! \return true if the bond satisfies the geometric properties for a double bond.

bool is_double_bond_geometry(int a_id, int b_id, hess::Molecule* mol) {
  Atom *a = mol->get_atom(a_id);
  Atom *b = mol->get_atom(b_id);
  if (a->valence > 3 || a->hybtype == 1 || b->valence > 3 || b->hybtype == 1) {
    return true;
  }
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_a_id = atom_ver.neiVertex(id);
    if (nei_a_id != b_id) {
      ver atom_ver_b = mol->getVertex(b_id);
      for (int id_b = atom_ver_b.neiBegin(); id_b != atom_ver_b.neiEnd(); id_b = atom_ver_b.neiNext(id_b)) {
        int nei_b_id = atom_ver_b.neiVertex(id_b);
        if (nei_b_id != a_id) {
          double torsion = fabs(get_torsion_angle(nei_a_id, a_id, b_id, nei_b_id, mol));
          if (torsion > 15.0 && torsion < 160.0) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

void assign_bond_orders(hess::Molecule* mol, bool after_aromatize) {
  vector <pair<int, double>> sorted_atoms;
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    Atom *a = mol->get_atom(a_id);
    double shortest_bond = 1.0e5;
    ver atom_ver = mol->getVertex(a_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      int bond_id = mol->findEdgeIndex(a_id, nei_id);
      Bond *bond = mol->get_bond(bond_id);
      if (mol->getAtomNumber(nei_id) != 1) {
        shortest_bond = min(shortest_bond, bond->length);
      }
    }
    double el_neg = GetElectroNeg(mol->getAtomNumber(a_id));
    pair<int, double> entry(a_id, el_neg * 1e6 + shortest_bond);
    sorted_atoms.push_back(entry);
  }
  sort(sorted_atoms.begin(), sorted_atoms.end(), sort_atom_z);
  for (int i = 0; i < sorted_atoms.size(); i++) {
    int a_id = sorted_atoms[i].first;
    Atom *a = mol->get_atom(a_id);
    if (((a->hybtype == 1) || mol->getVertex(a_id).degree() == 1)
            && a->valence + 2 <= GetMaxBonds(mol->getAtomNumber(a_id))) {
      if (has_bond_of_order(a_id, 2, mol) || has_bond_of_order(a_id, 3, mol)
              || (mol->getAtomNumber(a_id) == 7 && a->valence + 2 > 3)) {
        continue;
      }
      double max_el_neg = 0.0;
      double shortest_bond = 5000.0;
      int c_id = -1;
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        int bond_id = mol->findEdgeIndex(a_id, nei_id);
        Bond *bond = mol->get_bond(bond_id);
        Atom *b = mol->get_atom(nei_id);
        double curr_el_neg = GetElectroNeg(mol->getAtomNumber(nei_id));
        bool a_ring_property = false;
        bool b_ring_property = false;
        if (after_aromatize) {
          a_ring_property = a->in_ring;;
          b_ring_property = b->in_ring;;
        } else {
          a_ring_property = !a->in_ring;
          b_ring_property = !b->in_ring;
        }
        if (((b->hybtype == 1) || mol->getVertex(nei_id).degree() == 1)
                && b->valence + 2 <= GetMaxBonds(mol->getAtomNumber(nei_id))
                && (curr_el_neg > max_el_neg || (is_aprox(curr_el_neg, max_el_neg, 1.0e-6) && bond->length < shortest_bond))
                && a_ring_property && b_ring_property) {
          if (has_bond_of_order(nei_id, 2, mol) || has_bond_of_order(nei_id, 3, mol)
                  || (mol->getAtomNumber(nei_id) == 7 && b->valence + 2 > 3)) {
            continue;
          }
          double bond_length = bond->length;
          if (mol->getVertex(a_id).degree() == 1 || mol->getVertex(nei_id).degree() == 1) {
            double test_length = mol->get_correct_bond_rad(a_id) + mol->get_correct_bond_rad(nei_id);
            if (bond_length > 0.9 * test_length)
              continue;
          }
          shortest_bond = bond_length;
          max_el_neg = curr_el_neg;
          c_id = nei_id;
        }
      }
      if (c_id != -1) {
        bond_order_changed(a_id, c_id, 3, mol);
      }
    }
      // sp2-hybrid
    else if (((a->hybtype == 2) || mol->getVertex(a_id).degree() == 1)
            && a->valence + 1 <= GetMaxBonds(mol->getAtomNumber(a_id))
            ) {
      if (has_bond_of_order(a_id, 2, mol) || has_bond_of_order(a_id, 3, mol)
              || (mol->getAtomNumber(a_id) == 7 && a->valence + 1 > 3)) {
        continue;
      }

      if (a->in_ring && mol->getAtomNumber(a_id) == 16) {
        mol->charge_calculate();
        if (mol->total_charge > 1 && mol->getAtomCharge(a_id) == 0) {
          mol->setAtomCharge(a_id, 1);
        } else continue;
      }
      double max_el_neg = 0.0;
      double shortest_bond = 5000.0;
      int c_id = -1;
      Atom *c = nullptr;
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        int bond_id = mol->findEdgeIndex(a_id, nei_id);
        Bond *bond = mol->get_bond(bond_id);
        Atom *b = mol->get_atom(nei_id);
        double curr_el_neg = GetElectroNeg(mol->getAtomNumber(nei_id));
        bool a_ring_property = false;
        bool b_ring_property = false;
        if (after_aromatize) {
          a_ring_property = a->in_ring;
          b_ring_property = b->in_ring;
        } else {
          a_ring_property = !a->in_ring;
          b_ring_property = !b->in_ring;
        }
        if (((b->hybtype == 2) || mol->getVertex(nei_id).degree() == 1)
                && b->valence + 1 <= GetMaxBonds(mol->getAtomNumber(nei_id))
                && is_double_bond_geometry(a_id, nei_id, mol)
                && (curr_el_neg > max_el_neg || (is_aprox(curr_el_neg, max_el_neg, 1.0e-6)))
                && a_ring_property && b_ring_property) {
          if (has_bond_of_order(nei_id, 2, mol) || has_bond_of_order(nei_id, 3, mol)
                  || (mol->getAtomNumber(nei_id) == 7 && b->valence + 1 > 3)) {
            continue;
          }

          if (b->in_ring && mol->getAtomNumber(nei_id) == 16) {
            mol->charge_calculate();
            if (mol->total_charge > 1 && mol->getAtomCharge(nei_id) == 0) {
              mol->setAtomCharge(nei_id, 1);
            } else continue;
          }

          double bond_length = bond->length;
          if (mol->getVertex(a_id).degree() == 1 || mol->getVertex(nei_id).degree() == 1) {
            double test_length = mol->get_correct_bond_rad(a_id) + mol->get_correct_bond_rad(nei_id);
            if (bond_length > 0.93 * test_length)
              continue;
          }

          double difference = shortest_bond - bond_length;
          if ((difference > 0.1)
                  || ((difference > -0.01) &&
                  ((!a->in_ring || !c || !c->in_ring || b->in_ring)
                  || (a->in_ring && c && !c->in_ring && b->in_ring)))) {

            shortest_bond = bond_length;
            max_el_neg = curr_el_neg;
            c_id = nei_id;
            c = b;
          }
        }
      }
      if (c_id != -1) {
        bond_order_changed(a_id, c_id, 2, mol);
      }
    }
  }
}

bool has_bond_of_order(int a_id, int order, hess::Molecule* mol) {
  ver v = mol->getVertex(a_id);

  for (int id = v.neiBegin(); id != v.neiEnd(); id = v.neiNext(id))
    if (mol->getBondOrder(v.neiEdge(id)) == order)
      return true;

  return false;
}

void bond_order_changed(const int a_id, const int b_id, const int order, hess::Molecule* mol) {
  int bond_id = mol->findEdgeIndex(a_id, b_id);
  Atom *a = mol->get_atom(a_id);
  Atom *b = mol->get_atom(b_id);
  if (bond_id == -1)
    throw HessException("Attempt to change the order of a bond that doesn't exist");
  int bond_order = mol->getBondOrder(bond_id);

  if (order == 2 && bond_order == 2) {
    return;
  } else if (order == 3 && bond_order == 3) {
    return;
  } else if (order == 1 && bond_order == 1) {
    return;
  } else if (order == 2 && bond_order == 1) {
    b->valence++;
    a->valence++;
  } else if (order == 2 && bond_order == 3) {
    b->valence -= 2;
  } else if (order == 3 && bond_order == 1) {
    b->valence += 2;
    a->valence += 2;
  } else if (order == 3 && bond_order == 2) {
    b->valence++;
    a->valence++;
  } else if (order == 1 && bond_order == 2) {
    b->valence--;
    a->valence--;
  } else if (order == 1 && bond_order == 3) {
    b->valence -= 2;
    a->valence -= 2;
  }
  mol->setBondOrder(bond_id, order);
}

void delete_bonds_exceeding_valence(hess::Molecule* mol) {
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    Atom *a = mol->get_atom(a_id);
    int a_val = mol->getVertex(a_id).degree();
    int a_bonds_limit = GetMaxBonds(mol->getAtomNumber(a_id));
    while (a_val > a_bonds_limit) {
      fprintf(hessGetStream(), "Detected more bonds than are allowed!\n");
      double max_length = 0;
      int max_length_id = 0;
      ver atom_ver = mol->getVertex(a_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        int bond_id = mol->findEdgeIndex(a_id, nei_id);
        Bond* bond = mol->get_bond(bond_id);
        if (bond->length > max_length) {
          max_length = bond->length;
          max_length_id = bond_id;
        }
      }
      mol->remove_bond(max_length_id);
      a_val -= 1;
    }
  }
}

double distance(const hess::Atom & a, const hess::Atom & b) {
  return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

double distance_sqr(const hess::Atom & a, const hess::Atom & b) {
  return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
}

bool sort_z(const pair<Atom*, int>& a1, const pair<Atom*, int>& a2) {
  return a1.first->z < a2.first->z;
}

void determine_connectivity(hess::Molecule* mol) {
  int bonds_count = mol->edgeCount();
  if (bonds_count > 0) {
    mol->bonds.resize(mol->edgeEnd());
    memset(&(mol->bonds[0]), 0, sizeof (Bond) * mol->bonds.size());
  }
  set <pair<int, int>> ex_bonds;
  vector<pair < Atom*, int>> z_sorted_atoms;
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    Atom* a = mol->get_atom(a_id);
    z_sorted_atoms.push_back(make_pair(a, a_id));
  }
  double max_rad = 0;
  double up_border = 0.45;
  double down_border = 0.16;
  sort(z_sorted_atoms.begin(), z_sorted_atoms.end(), sort_z);
  vector<double> rads;
  for (int i = 0; i < z_sorted_atoms.size(); i++) {
    double rad = GetCovalentRad(mol->getAtomNumber(z_sorted_atoms[i].second));
    max_rad = max(rad, max_rad);
    rads.push_back(rad);
  }
  for (int i = 0; i < z_sorted_atoms.size(); i++) {
    Atom* a = z_sorted_atoms[i].first;
    int a_id = z_sorted_atoms[i].second;
    int b_id;
    double maxcutoff = (rads[i] + max_rad + up_border) * (rads[i] + max_rad + up_border);
    Atom* b = nullptr;
    for (int j = i + 1; j < z_sorted_atoms.size(); j++) {
      b = z_sorted_atoms[j].first;
      b_id = z_sorted_atoms[j].second;
      if ((a_id == b_id) || (ex_bonds.find(make_pair(b_id, a_id)) != ex_bonds.end())) {
        continue;
      }
      if (a->valence >= GetMaxBonds(mol->getAtomNumber(a_id)) ||
              b->valence >= GetMaxBonds(mol->getAtomNumber(b_id))) {
        continue;
      }
      if (mol->getAtomNumber(a_id) == 1 && mol->getAtomNumber(b_id) == 1) {
        continue;
      }
      double cutoff = (rads[i] + rads[j] + up_border) * (rads[i] + rads[j] + up_border);
      double zd = (a->z - b->z) * (a->z - b->z);
      if (zd > maxcutoff)
        break;
      double xd = (a->x - b->x) * (a->x - b->x);
      if (xd > cutoff)
        continue;
      xd += (a->y - b->y) * (a->y - b->y);
      if (xd > cutoff)
        continue;
      xd += zd;
      if (xd > cutoff || xd < down_border)
        continue;
      double length = distance(*a, *b);
      a->valence++;
      b->valence++;
      int bond_id = mol->findEdgeIndex(a_id, b_id);
      if (bond_id >= 0)
        mol->bonds[bond_id].length = length;
      else
        mol->add_bond(a_id, b_id, 1, length);
      ex_bonds.insert(make_pair(a_id, b_id));
    }
  }
}
