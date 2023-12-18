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
 * https://github.com/mwojcikowski/smina/blob/master/src/lib/everything.h
 * Copyright (c) 2006-2010, The Scripps Research Institute
 ***********************************************************************/

#include "features.h"
#include "elements.h"

using namespace hess;

// Convenience function used by the DeleteHydrogens methods

bool is_suppressible_hydrogen(hess::Molecule* mol, int atom_id) {
  if (mol->getAtomIsotope(atom_id) == 0 && nonhydrobond_count(mol, atom_id) == 1 && mol->getAtomCharge(atom_id) == 0)
    return true;
  else
    return false;
}

bool is_rot_bond(hess::Molecule *mol, int bond_id) {
  Bond* bond = mol->get_bond(bond_id);
  if (!(mol->getBondOrder(bond_id) == 1) || bond_is_amide(mol, bond_id) || bond->in_ring) {
    return false;
  }
  if ((nonhydrobond_count(mol, mol->bond_atom_beg(bond_id)) == 1) || (nonhydrobond_count(mol, mol->bond_atom_end(bond_id)) == 1)) {
    return false;
  }
  return true;
}

bool atom_is_nitro_oxygen(hess::Molecule *mol, int a_id) {
  //atom is connected to a nitrogen that has a total
  //of 2 attached free oxygens
  if (mol->getAtomNumber(a_id) != 8)
    return false;
  if (nonhydrobond_count(mol, a_id) != 1) {
    return false;
  }
  int atom_id = -1;
  ver atom_ver = mol->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (mol->getAtomNumber(nei_id) == 7) {
      atom_id = nei_id;
      break;
    }
  }
  if (atom_id != -1)
    return false;
  if (count_free_oxygens(mol, atom_id) != 2)
    return false;
  return true;
}

bool atom_is_sulfone_oxygen(hess::Molecule *mol, int atom_id) {
  if (mol->getAtomNumber(atom_id) != 8)
    return false;
  if (nonhydrobond_count(mol, atom_id) != 1) {
    return false;
  }
  int nbr_id = -1;
  ver atom_ver = mol->getVertex(atom_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (mol->getAtomNumber(nei_id) == 16) {
      nbr_id = nei_id;
      break;
    }
  }
  if (nbr_id == -1) {
    return false;
  }
  // check for sulfate
  if (count_free_oxygens(mol, nbr_id) != 2)
    return false;

  ver atom_ver2 = mol->getVertex(nbr_id);
  for (int id2 = atom_ver2.neiBegin(); id2 != atom_ver2.neiEnd(); id2 = atom_ver2.neiNext(id2)) {
    int nei2_id = atom_ver2.neiVertex(id2);
    Atom *nei2 = mol->get_atom(nei2_id);
    if (mol->getAtomNumber(nei2_id) == 7) {
      return false;
    }
  }
  return true; // true sulfone
}

bool bond_is_carbonyl(hess::Molecule *mol, int bgn_id, int end_id) {
  const int bond_id = mol->findEdgeIndex(bgn_id, end_id);
  const Bond *bond = mol->get_bond(bond_id);
  if (mol->getBondOrder(bond_id) != 2)
    return false;

  if ((mol->getAtomNumber(bgn_id) == 6 && mol->getAtomNumber(end_id) == 8) ||
          (mol->getAtomNumber(bgn_id) == 8 && mol->getAtomNumber(end_id) == 6))
    return true;

  return false;
}


//Is the bond an ester link (i.e., between a carbonyl C and an O)?

bool bond_is_ester(hess::Molecule *mol, int bgn_id, int end_id) {
  int a1 = -1, a2 = -1;
  const int bond_id = mol->findEdgeIndex(bgn_id, end_id);

  if (mol->getAtomNumber(bgn_id) == 6 && mol->getAtomNumber(end_id) == 8) {
    a1 = bgn_id;
    a2 = end_id;
  }

  if (mol->getAtomNumber(bgn_id) == 8 && mol->getAtomNumber(end_id) == 6) {
    a1 = end_id;
    a2 = bgn_id;
  }

  if (a1 == -1 || a2 == -1)
    return false;
  if (mol->getBondOrder(bond_id) != 1)
    return false;
  ver atom_ver = mol->getVertex(a1);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (bond_is_carbonyl(mol, a1, nei_id)) {
      return true;
    }
  }
  return false;
}

bool bond_is_amide(hess::Molecule *mol, int bond_id) {
  int c_id = -1, n_id = -1;
  Atom *n = nullptr;
  Bond* bond = mol->get_bond(bond_id);
  int bgn_id = mol->bond_atom_beg(bond_id);
  int end_id = mol->bond_atom_end(bond_id);
  Atom *bgn = mol->get_atom(bgn_id);
  Atom *end = mol->get_atom(end_id);

  // Look for C-N bond
  if (mol->getAtomNumber(bgn_id) == 6 && mol->getAtomNumber(end_id) == 7) {
    c_id = bgn_id;
    n_id = end_id;
    n = end;
  } else if (mol->getAtomNumber(bgn_id) == 7 && mol->getAtomNumber(end_id) == 6) {
    c_id = end_id;
    n_id = bgn_id;
    n = bgn;
  }
  if (c_id == -1 || n_id == -1) return (false);
  if (mol->getBondOrder(bond_id) != 1) return (false);
  if (n && mol->getVertex(n_id).degree() != 3) return false; // must be a degree 3 nitrogen

  // Make sure C is attached to =O
  ver atom_ver = mol->getVertex(c_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    if (bond_is_carbonyl(mol, c_id, nei_id))
      return true;
  }

  return false;
}



// IsHbondAcceptor
// This function has a finer grain than the original
// implementation, checking also the neighbor atoms.
// Accordingly to these rules, the function will return:
//
//    aliph-O-aliph ether   -> true   [1]
//    hydroxy O-sp3         -> true   [1]
//    aro-O-aliph ether     -> true   [1]
//    ester O-sp2           -> true   [1]
//    sulfate O (R-SO3)     -> true   [2]
//    sulfoxyde O (R-SO-R)  -> true   [2]
//    organoboron-F (R-BF3) -> true   [3]
//    ester O-sp3           -> false  [1]
//    sulfone (R1-SO2-R2 )  -> false  [2]
//    aro-O-aro             -> false  [1]
//    aromatic O            -> false  [1]
//    O-nitro               -> false  [2]
//    organic F (R-F)       -> false  [4]
//

bool atom_is_Hbond_acceptor(hess::Molecule *mol, int atom_id) {
  Atom *atom = mol->get_atom(atom_id);
  if (mol->getAtomNumber(atom_id) == 8) {
    unsigned int aro_count = 0;

    if (atom_is_nitro_oxygen(mol, atom_id)) { // maybe could be a bool option in the function?
      return false;
    }
    if (atom->aromatic) { // aromatic oxygen (furan) (NO)
      return false;
    }
    if (atom_is_sulfone_oxygen(mol, atom_id)) { // sulfone (NO)
      return false;
    }
    ver atom_ver = mol->getVertex(atom_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      Atom *nbr = mol->get_atom(nei_id);
      if (nbr->aromatic) {
        aro_count += 1;
        if (aro_count == 2) { // aromatic ether (aro-O-aro) (NO)
          return false;
        }
      } else {
        if (mol->getAtomNumber(nei_id) == 1) // hydroxyl (YES)
          return true;
        if ((bond_is_ester(mol, nei_id, atom_id)) && (!(is_carboxyl_oxygen(mol, atom_id))))
          return false;
      }
    }
    return true; // any other oxygen
  } // oxygen END

  // fluorine
  if (mol->getAtomNumber(atom_id) == 9) {
    // organic fluorine (NO)
    ver atom_ver = mol->getVertex(atom_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      Atom *nbr = mol->get_atom(nei_id);
      if (mol->getAtomNumber(nei_id) == 6)
        return false;
      else
        return true;
    }
  }
  if (mol->getAtomNumber(atom_id) == 7) {
    if (!((mol->getVertex(atom_id).degree() == 4 && atom->hybtype == 3)
            || (mol->getVertex(atom_id).degree() == 3 && atom->hybtype == 2)))
      return true;
  }
  if (mol->getAtomNumber(atom_id) == 16 && mol->getAtomCharge(atom_id) == -1) {
    return true;
  }
  return false;
}

bool bonded_to_HD(hess::Molecule *mol, const int atom_id) {
  ver atom_ver = mol->getVertex(atom_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    Atom *nbr = mol->get_atom(nei_id);
    if (nbr->atom_type == specific_atom_type::PolarHydrogen) {
      return true;
    }
  }
  return false;
}

bool bonded_to_heteroatom(hess::Molecule *mol, const int atom_id) {
  ver atom_ver = mol->getVertex(atom_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    Atom *nbr = mol->get_atom(nei_id);
    if (is_heteroatom(nbr->atom_type)) {
      return true;
    }
  }
  return false;
}

smt adjust_atom_type(smt t, bool Hbonded, bool heteroBonded) {
  using namespace specific_atom_type;
  switch (t) {
    case AliphaticCarbonXSHydrophobe: // C_C_C_H, //hydrophobic according to xscale
    case AliphaticCarbonXSNonHydrophobe: //C_C_C_P,
      return heteroBonded ? AliphaticCarbonXSNonHydrophobe : AliphaticCarbonXSHydrophobe;
    case AromaticCarbonXSHydrophobe: //C_A_C_H,
    case AromaticCarbonXSNonHydrophobe: //C_A_C_P,
      return heteroBonded ? AromaticCarbonXSNonHydrophobe : AromaticCarbonXSHydrophobe;
    case NitrogenXSDonor: //N_N_N_D,
    case specific_atom_type::Nitrogen: //N_N_N_P, no hydrogen bonding
      return Hbonded ? NitrogenXSDonor : specific_atom_type::Nitrogen;
    case NitrogenXSDonorAcceptor: //N_NA_N_DA, also an autodock acceptor
    case NitrogenXSAcceptor: //N_NA_N_A, also considered an acceptor by autodock
      return Hbonded ? NitrogenXSDonorAcceptor : NitrogenXSAcceptor;
    case OxygenXSDonor: //O_O_O_D,
    case specific_atom_type::Oxygen: //O_O_O_P,
      return Hbonded ? specific_atom_type::OxygenXSDonor : specific_atom_type::Oxygen;
    case OxygenXSDonorAcceptor: //O_OA_O_DA, also an autodock acceptor
    case OxygenXSAcceptor: //O_OA_O_A, also an autodock acceptor
      return Hbonded ? OxygenXSDonorAcceptor : OxygenXSAcceptor;
    default:
      return t;
  }

}

void adjust_atom_types(hess::Molecule *mol) {
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    Atom* atom = mol->get_atom(a_id);
    atom->atom_type = adjust_atom_type(atom->atom_type, bonded_to_HD(mol, a_id), bonded_to_heteroatom(mol, a_id));
  }
}

void assign_atom_types(hess::Molecule* mol, const string& scoring_type_str) {
  for (size_t i = 0u; i < specific_atom_type::NumTypes; ++i)
    if (scoring_type_str == "vinardo")
      specific_atom_type::data[i] = specific_atom_type::vinardo_data[i];
    else
      specific_atom_type::data[i] = specific_atom_type::default_data[i];
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    Atom* atom = mol->get_atom(a_id);
    int atomic_num = mol->getAtomNumber(a_id);
    if (atomic_num == 0)
      throw HessException("unexpected zero atomic number");
    const char *element_name = GetSymbol(atomic_num);
    char element_name_final[3];
    element_name_final[2] = '\0';
    if (atomic_num == 1) {
      element_name_final[0] = 'H';
      element_name_final[1] = 'D';
    } else if ((atomic_num == 6) && (atom->aromatic)) {
      element_name_final[0] = 'A';
      element_name_final[1] = '\0';
    } else if (atomic_num == 8) {
      element_name_final[0] = 'O';
      element_name_final[1] = 'A';
    } else if ((atomic_num == 7) && (atom_is_Hbond_acceptor(mol, a_id))) {
      element_name_final[0] = 'N';
      element_name_final[1] = 'A';
    } else if ((atomic_num == 16) && (atom_is_Hbond_acceptor(mol, a_id))) {
      element_name_final[0] = 'S';
      element_name_final[1] = 'A';
    } else {
      if (!isalnum(element_name[0])) {
        element_name_final[0] = '\0';
      } else {
        element_name_final[0] = element_name[0];
      }
      if (!isalnum(element_name[1])) {
        element_name_final[1] = '\0'; //null terminate
      } else {
        element_name_final[1] = element_name[1];
      }
    }
    smt sm = string_to_atom_type(element_name_final);
    if (sm >= specific_atom_type::NumTypes) {
      int g = 1;
    }
    assert(sm < specific_atom_type::NumTypes);
    atom->atom_type = sm;
  }
}

double gaussian(double x, double width) {
  return exp(-sqr(x / width));
}

double smooth_div(double x, double y) {
  if (std::abs(x) < epsilon_fl) return 0;
  if (std::abs(y) < epsilon_fl) return ((x * y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
  return x / y;
}

double slope_step(double x_bad, double x_good, double x) {
  if (x_bad < x_good) {
    if (x <= x_bad)
      return 0;
    if (x >= x_good)
      return 1;
  } else {
    if (x >= x_bad)
      return 0;
    if (x <= x_good)
      return 1;
  }
  return (x - x_bad) / (x_good - x_bad);
}

double slope_step_deriv(double x_bad, double x_good, double x) {
  if (x_bad < x_good) {
    if (x <= x_bad)
      return 0;
    if (x >= x_good)
      return 0;
  } else {
    if (x >= x_bad)
      return 0;
    if (x <= x_good)
      return 0;
  }
  return 1.0 / (x_good - x_bad);
}

double optimal_distance(smt xs_t1, smt xs_t2) {
  return xs_radius(xs_t1) + xs_radius(xs_t2);
}

ConfIndependentInputs::ConfIndependentInputs() :
num_tors(0), num_rotors(0), num_heavy_atoms(0),
num_hydrophobic_atoms(0), ligand_max_num_h_bonds(0), num_ligands(
0),
ligand_lengths_sum(0) {
}

ConfIndependentInputs::ConfIndependentInputs(hess::Molecule* ligand) {
  num_tors = 0;
  num_rotors = 0;
  num_heavy_atoms = 0;
  num_hydrophobic_atoms = 0;
  ligand_max_num_h_bonds = 0;
  num_ligands = 1;
  ligand_lengths_sum = 0;
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    hess::Atom* a = ligand->get_atom(a_id);
    if (!is_hydrogen(a->atom_type)) {
      unsigned ar = atom_rotors(ligand, a_id);

      num_tors += 0.5 * ar;

      if (ar > 2)
        num_rotors += 0.5;
      else
        num_rotors += 0.5 * ar;

      ++num_heavy_atoms;
      if (xs_is_hydrophobic(a->atom_type))
        ++num_hydrophobic_atoms;

      if (xs_is_acceptor(a->atom_type) || xs_is_donor(a->atom_type))
        ++ligand_max_num_h_bonds;
    }
  }
}


// the number of rotatable bonds to heavy ligand atoms

unsigned ConfIndependentInputs::atom_rotors(hess::Molecule* ligand, const int a_id) {
  unsigned acc = 0;
  ver atom_ver = ligand->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    Atom *b = ligand->get_atom(nei_id);
    Bond *bond = ligand->get_bond(ligand->findEdgeIndex(a_id, nei_id));
    if (!is_hydrogen(b->atom_type) && bond->rotatable && this->num_bonded_heavy_atoms(ligand, nei_id) > 1) {
      if (this->num_bonded_heavy_atoms(ligand, a_id) > 1) {
        ++acc;
      }
    }
  }
  return acc;
}

unsigned ConfIndependentInputs::num_bonded_heavy_atoms(hess::Molecule* ligand, int a_id) {
  unsigned acc = 0;
  ver atom_ver = ligand->getVertex(a_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    Atom *b = ligand->get_atom(nei_id);
    if (!is_hydrogen(b->atom_type)) {
      ++acc;
    }
  }
  return acc;
}