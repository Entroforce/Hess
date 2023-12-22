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

/* Some contents of this file was borrowed from  Openbabel repository
 * https://github.com/openbabel/openbabel/blob/master/data/bondtyp.txt
 * Copyright (c) 2002-2005 by Geoffrey R Hutchison                           
 * Part of the Open Babel package, under the GNU General Public License (GPL)
 */

#include "bond_order_determination.h"
#include "molecule/elements.h"

#include "base_cpp/scanner.h"
#include "molecule/smiles_loader.h"
#include "molecule/query_molecule.h"
#include "molecule/molecule_substructure_matcher.h"


const vector<vector<string>> func_groups = {
  {"[X2,X3]1[#6]([#7D3]2)[#6][#6][#6]2[X2,X3][#6]([#7D3]3)[#6][#6][#6]3[X2,X3][#6]([#7D3]4)[#6][#6][#6]4[X2,X3][#6]([#7D3]5)[#6][#6][#6]51",
    "0 1 2 1 2 1 1 3 1 3 4 2 4 5 1 5 2 1 5 6 2 6 7 1 7 8 2 7 9 1 9 10 2 10 11 1 11 8 1 11 12 2 12 13 1 13 14 1 13 15 2 15 16 1 16 17 2 17 14 1 17 18 1 18 19 2 19 20 1 19 21 1 21 22 2 22 23 1 23 20 2",
    "!"},
  {"[X2,X3]1[#6]([#7D3]2)[#6][#6][#6]2[X2,X3][#6]([#7]3)[#6][#6][#6]3[X2,X3][#6]([#7D3]4)[#6][#6][#6]4[X2,X3][#6]([#7]5)[#6][#6][#6]51",
    "0 1 2 1 2 1 1 3 1 3 4 2 4 5 1 5 2 1 5 6 2 6 7 1 7 8 2 7 9 1 9 10 2 10 11 1 11 8 1 11 12 2 12 13 1 13 14 1 13 15 2 15 16 1 16 17 2 17 14 1 17 18 1 18 19 2 19 20 1 19 21 1 21 22 2 22 23 1 23 20 2",
    "!"},
  {"[X2,X3]1[#6]([#7]2)[#6][#6][#6]2[X2,X3][#6]([#7]3)[#6][#6][#6]3[X2,X3][#6]([#7]4)[#6][#6][#6]4[X2,X3][#6]([#7]5)[#6][#6][#6]51",
    "0 1 2 1 2 1 1 3 1 3 4 2 4 5 1 5 2 1 5 6 2 6 7 1 7 8 2 7 9 1 9 10 2 10 11 1 11 8 1 11 12 2 12 13 1 13 14 1 13 15 2 15 16 1 16 17 2 17 14 1 17 18 1 18 19 2 19 20 1 19 21 1 21 22 2 22 23 1 23 20 2",
    "!"},
  {"[#7D2][#7D2][#7D1]", "0 1 2 1 2 2", "_1_"}, //# Azide
  //    {"[#8D1][#7D3]([#8D1])*", "0 1 2 1 2 2 1 3 1", "_2_"}, //# Nitro
  {"[#16D4]([#8D1])([#8D1])([*!#8])([*!#8])", "0 1 2 0 2 2 0 3 1 0 4 1", "___"}, //# Sulfones
  {"[#16D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])", "0 1 2 0 2 2 0 3 1 0 4 1", "___"}, //# Sulfates
  {"[#16D4]([#16D1])([#8D1])([#8-,#8])([#8-,#8])", "0 1 2 0 2 2 0 3 1 0 4 1", "___"}, //# Thiosulfates
  {"[#16D3]([#8D1])([*!#8])([*!#8])", "0 1 2 0 2 1 0 3 1", "__"}, //# Sulfoxides
  {"[#16D3]([#8D1])([#8D1-])([#8D1-])", "0 1 2 0 2 1 0 3 1", "__"}, //# Sulfite
  {"[#16D3]([#8D1])([#8D1])([#8D1])", "0 1 2 0 2 2 0 3 2", "____"}, //# Sulfur trioxide
  {"[#16D3]([#8D1])([#8])([#8])", "0 1 2 0 2 1 0 3 1", "__"}, //# Sulfites
  {"[#16D2]([#8D1])([#16D1])", "0 1 2 0 2 2", "___"}, //# Disulfur monoxide
  {"[#16D2]([#8D1])([*!#8])", "0 1 2 0 2 1", "__"}, //# Sulfmonoxides
  {"[#16D2]([#8D1])([#8D1])", "0 1 2 0 2 2", "___"}, //# Sulfur dioxide
  {"[#15D3]([#8D1])([#8D1])([#8D2])", "0 1 2 0 2 2 0 3 1", "__"}, //#Phosphite
  //    {"[#7D2]([#8D1])([#1])", "0 1 2 0 2 1", "__"},
  {"[#15D4]([#8D1])(*)(*)(*)", "0 1 2 0 2 1 0 3 1 0 4 1", "__"}, //# Phosphone
  {"[#6D3]([#8D1])([#8])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# Carboxylic Acid, ester, etc.
  {"[#8D1][#6D2][#8D1]", "0 1 2 1 2 2", "_1_"}, //# Carbon dioxide
  {"[#6D3]([#8D1;!-])([#7])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# Amide C(=O)N  - no negative charge on O (2aix_u1k.sdf)
  {"[#34D3]([#8D1])([#8])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# Seleninic acid Se(=O)OH
  {"[#6D3]([#8D1])([#16])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# Thioacid / Thioester C(=O)S
  {"[#6D3]([#16D1])([#16])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# dithioacid / dithioester C(=S)S
  {"[#6D3]([#16D1])([#7])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# thioamide C(=S)N # avoid aromatics (pdb_ligands_sdf/1yry_msg.sdf)
  {"[#6][#6D2][#6]", "0 1 2 1 2 2", "212"}, //# allene C=C=C
  {"[#6][#6D2][#8D1]", "0 1 2 1 2 2", "21_"}, //# ene-one C=C=O
  {"[#6D1][#7D2]*", "0 1 3 1 2 1", "_1"}, //# isonitrile / isocyano
  {"[Nv2R][#6v3][#8v2]", "0 1 2 1 2 1", "_2_"}, //# NR2 in ring with hybridized carbon neighbors, do not apply to aromatic N
  {"[Nv2R][#6v3][Nv2]", "0 1 3 1 2 1", "_2_"}, //# valence insted of degree used to fix pdb_ligands_sdf/3dcv_55e.sdf
  {"[#6D3;!R]([#7D2;!R])([#7D1;!R])[#7D1;!R]", "0 1 2 0 2 1 0 3 1", "2____"}, //# if three N are present in R-N-guanidine-ish
  {"[#6D3;!R]([#7D1H0;!R])([#7;!R])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# guanidinium and amidine -C(=NH)NH2 without hydrogens
  {"[#6D3;!R]([#7D2H1;!R])([#7;!R])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# (this can normally be figured out, but is needed to avoid matching the next SMARTS)
  {"[#6D3;!R]([#7D3H2;!R])([#7;!R])*", "0 1 2 0 2 1 0 3 1", "2__"}, //# and also with more hydrogens than normal (protonated)
  {"[#6D3;!R]([#1,#6])([#1,#6])[#7D3;!R]([#1])[#6]", "0 1 1 0 2 1 0 3 2 3 4 1 3 5 1", "2__2__"}, //# Schiff base, protonated
};

using namespace hess;
using namespace indigo;

class BondTriplet {
public:

  BondTriplet(int fir, int sec, int _order) : first(fir), second(sec), order(_order) {
  }
  int first;
  int second;
  int order;
};

class FGroup {
public:

  FGroup(vector<BondTriplet>&& bonds, const string& smarts_format, vector <int>&& hybs) {
    this->bonds = bonds;
    this->smarts_format = smarts_format;
    this->hybrs = hybs;
  }

  FGroup(FGroup&& other) {
    this->bonds = move(other.bonds);
    this->smarts_format = move(other.smarts_format);
    this->hybrs = move(other.hybrs);
    other.bonds.clear();
    other.smarts_format = "";
    other.hybrs.clear();
  }

  FGroup(const FGroup& other) {
    this->bonds = other.bonds;
    this->smarts_format = other.smarts_format;
    this->hybrs = other.hybrs;
  }

  string smarts_format;
  vector<BondTriplet> bonds;
  vector <int> hybrs;
};

static void _parse(vector <FGroup> & groups) {
  for (const vector<string> &func_group : func_groups) {
    if (func_group.size() != 3) {
      throw HessException("Incorrect representation of functional group patterns");
    }
    vector <BondTriplet> bonds;
    string atoms_buff = func_group[0];
    string bonds_buff = func_group[1];
    string sp_buff = func_group[2];
    int index = 0;
    string scope = "";
    bool open_round = false;
    while (index != atoms_buff.size()) {
      if (atoms_buff[index] == '(') {
        if (open_round) {
          throw HessException("Groups parser error. Detected open without closing");
        }
        open_round = true;
        scope += atoms_buff[index];
        index += 1;
        continue;
      } else if (atoms_buff[index] == ']' && open_round) {
        scope += atoms_buff[index];
        index += 1;
        continue;
      }
      else if (atoms_buff[index] == ']' && !open_round) {
        scope += atoms_buff[index];
        index += 1;
        scope = "";
        continue;
      } else if (atoms_buff[index] == '*' && !open_round) {
        scope += atoms_buff[index];
        index += 1;
        scope = "";
        continue;
      } else if (atoms_buff[index] == ')' && open_round) {
        scope += atoms_buff[index];
        index += 1;
        scope = "";
        open_round = false;
        continue;
      } else if (atoms_buff[index] == ')' && !open_round) {
        throw HessException("Groups parser error. Detected close without opening");
      }
      scope += atoms_buff[index];
      index += 1;
    }
    index = 0;
    while (index != bonds_buff.size()) {
      string bond_buff = "";
      int bonds_buf_len = bonds_buff.size();
      while (bonds_buff[index] != ' ' && index < bonds_buf_len) {
        bond_buff += bonds_buff[index];
        index += 1;
      }
      int first_index = atoi(bond_buff.c_str());
      index += 1;

      bond_buff = "";
      while (bonds_buff[index] != ' ' && index < bonds_buf_len) {
        bond_buff += bonds_buff[index];
        index += 1;
      }
      int second_index = atoi(bond_buff.c_str());
      index += 1;

      bond_buff = "";
      while (bonds_buff[index] != ' ' && index != bonds_buff.size()) {
        bond_buff += bonds_buff[index];
        index += 1;
      }
      int order = atoi(bond_buff.c_str());
      if (index != bonds_buff.size()) index += 1;

      BondTriplet triplet(first_index, second_index, order);
      bonds.push_back(triplet);
    }
    index = 0;
    vector <int> hybrs;
    if (sp_buff[0] == '!') {
      hybrs = {};
    } else {
      while (index != sp_buff.size()) {
        if (sp_buff[index] == '_') {
          hybrs.push_back(0);
        } else {
          char order = sp_buff[index];
          if (order == '1') hybrs.push_back(1);
          else if (order == '2') hybrs.push_back(2);
          else if (order == '3') hybrs.push_back(3);
        }
        index += 1;
      }
    }
    groups.push_back(FGroup(move(bonds), atoms_buff, move(hybrs)));
  }
}

void assign_functional_group_bonds(const vector <FGroup>& groups, hess::Molecule* mol) {
  for (int index = 0; index < groups.size(); index++) {
    const FGroup& group = groups[index];
    const vector<BondTriplet> &bonds = group.bonds;
    const vector<int> &hybrs = group.hybrs;

    indigo::BufferScanner scanner(group.smarts_format.c_str());
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int i;

      for (i = 0; i < hybrs.size(); i++)
        if (mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[i]])->hybtype != hybrs[i])
          break;
      if (i < hybrs.size())
        continue;

      for (const BondTriplet& triplet : bonds) {
        int a_id = molcopy_mapping[matcher.getQueryMapping()[triplet.first]];
        int b_id = molcopy_mapping[matcher.getQueryMapping()[triplet.second]];
        bond_order_changed(a_id, b_id, triplet.order, mol);
      }
    }
  }
}

void functional_groups(hess::Molecule* mol) {
  vector <FGroup> groups;
  _parse(groups);
  assign_functional_group_bonds(groups, mol);

  {
    indigo::BufferScanner scanner("[#8D1][#6](*)(*)");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      hess::Atom *a = mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[0]]);
      hess::Atom *b = mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[1]]);
      hess::Atom *c1 = mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[2]]);
      hess::Atom *c2 = mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[3]]);

      double ang1 = get_angle(*a, *b, *c1);
      double ang2 = get_angle(*a, *b, *c2);
      double av_angle = (ang1 + ang2) / 2;
      double dist = distance(*a, *b);
      if (av_angle > 115 && av_angle < 150 && dist < 1.28) {
        if (!has_bond_of_order(matcher.getQueryMapping()[0], 2, mol)) {
          bond_order_changed(matcher.getQueryMapping()[0], matcher.getQueryMapping()[1], 2, mol);
        }
      }
    }
  }

  {
    indigo::BufferScanner scanner("[#8D1][#7D3]([#8D1])*");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int res_1 = molcopy_mapping[matcher.getQueryMapping()[0]];
      int res_2 = molcopy_mapping[matcher.getQueryMapping()[1]];
      int res_3 = molcopy_mapping[matcher.getQueryMapping()[2]];

      if (mol->get_atom(res_2)->hybtype != 2) {
        continue;
      }

      bond_order_changed(res_1, res_2, 2, mol);
      mol->setAtomCharge(res_2, 1);
      mol->setAtomCharge(res_3, -1);
    }
  }


  {
    indigo::BufferScanner scanner("[#16D1][#6](*)(*)");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int res_1 = molcopy_mapping[matcher.getQueryMapping()[0]];
      int res_2 = molcopy_mapping[matcher.getQueryMapping()[1]];
      int res_3 = molcopy_mapping[matcher.getQueryMapping()[2]];
      int res_4 = molcopy_mapping[matcher.getQueryMapping()[3]];

      hess::Atom *a = mol->get_atom(res_1);
      hess::Atom *b = mol->get_atom(res_2);
      hess::Atom *c1 = mol->get_atom(res_3);
      hess::Atom *c2 = mol->get_atom(res_4);
      double ang1 = get_angle(*a, *b, *c1);
      double ang2 = get_angle(*a, *b, *c2);
      double av_angle = (ang1 + ang2) / 2;
      double dist = distance(*a, *b);
      if (av_angle > 115 && av_angle < 150 && dist < 1.72) {
        if (!has_bond_of_order(res_1, 2, mol)) {
          bond_order_changed(res_1, res_2, 2, mol);
        }
      }
    }
  }


  {
    indigo::BufferScanner scanner("[#8,#16;D1][#6D2][#7D2]");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int res_1 = molcopy_mapping[matcher.getQueryMapping()[0]];
      int res_2 = molcopy_mapping[matcher.getQueryMapping()[1]];
      int res_3 = molcopy_mapping[matcher.getQueryMapping()[2]];

      hess::Atom *a = mol->get_atom(res_1);
      hess::Atom *b = mol->get_atom(res_2);
      hess::Atom *c1 = mol->get_atom(res_3);

      double ang1 = get_angle(*a, *b, *c1);
      double dist1 = distance(*a, *b);
      double dist2 = distance(*b, *c1);
      bool distok = false;
      if (mol->getAtomNumber(res_1) == indigo::ELEM_O) {
        distok = dist1 < 1.28;
      } else {
        distok = dist1 < 1.72;
      }

      if (ang1 > 150 && distok && dist2 < 1.34) {
        bond_order_changed(res_1, res_2, 2, mol);
        bond_order_changed(res_3, res_2, 2, mol);
      }
    }
  }

  {
    indigo::BufferScanner scanner("[#6D3][#7D2][#8D2]");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int res_1 = molcopy_mapping[matcher.getQueryMapping()[0]];
      int res_2 = molcopy_mapping[matcher.getQueryMapping()[1]];
      int res_3 = molcopy_mapping[matcher.getQueryMapping()[2]];

      hess::Atom *a = mol->get_atom(res_1);
      hess::Atom *b = mol->get_atom(res_2);
      hess::Atom *c1 = mol->get_atom(res_3);

      double ang1 = get_angle(*a, *b, *c1);
      double dist = distance(*a, *b);
      if (ang1 > 110 && ang1 < 150 && dist < 1.4) {
        if (!has_bond_of_order(res_1, 2, mol)) {
          bond_order_changed(res_1, res_2, 2, mol);
        }
      }
    }
  }

  {
    indigo::BufferScanner scanner("[#8D1][#7D3r6]");
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int res_1 = molcopy_mapping[matcher.getQueryMapping()[0]];
      int res_2 = molcopy_mapping[matcher.getQueryMapping()[1]];

      hess::Atom *a = mol->get_atom(res_1);
      hess::Atom *b = mol->get_atom(res_2);

      double av_angle = 0.0;
      vector<hess::Atom*> neis;
      ver atom_ver = mol->getVertex(res_2);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        if (nei_id != res_1) {
          hess::Atom *nei = mol->get_atom(nei_id);
          neis.push_back(nei);
          double ang = get_angle(*a, *b, *nei);
          av_angle += ang;
        }
      }
      assert(neis.size() == 2);
      double ang = get_angle(*neis[0], *b, *neis[1]);
      av_angle += ang;
      av_angle /= 3;
      double dist = distance(*a, *b);
      if (av_angle > 110 && av_angle < 150 && dist < 1.35) {
        mol->setAtomCharge(res_1, -1);
        mol->setAtomCharge(res_2, 1);
      }
    }
  }
}