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

#include "base_cpp/scanner.h"
#include "molecule/smiles_loader.h"
#include "molecule/query_molecule.h"
#include "base_cpp/array.h"
#include "molecule/molecule_substructure_matcher.h"

using namespace std;
using namespace hess;

const vector<vector<string>> hyb_patterns = {
  {"[D4]", "3"},
  {"[D5]", "5"},
  {"[D6]", "6"},
  {"[C]", "3"},
  {"[c,$(C=*)]", "2"},
  {"[$(C#*),$(C(=*)=*)]", "1"},
  {"[N]", "3"},
  //    {"[ND1,ND2,ND3]a", "2"},
  {"[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]", "2"},
  {"[$(N#*),$([ND2](=*)=*)]", "1"},
  {"[O]", "3"},
  {"[o,$(O=*),$(O[#6,#7,#8]=,:*)]", "2"},
  {"[$([#8D1][#6][#8D1])]", "2"},
  {"[$(O#*)]", "1"},
  {"[P]", "3"},
  {"[#15;$([PD1]=*)]", "2"},
  {"[PD5]", "5"},
  {"[Pv5]", "5"},
  {"[S]", "3"},
  {"[#16;s,$([SD1]=*)]", "2"},
  {"[SD6]", "6"},
  {"[B]", "2"},
  {"[BD4]", "3"},
  {"[Al]", "2"},
  {"[Ga]", "2"},
  {"[In]", "2"},
  {"[Tl]", "2"},
  {"[AlD4]", "3"},
  {"[Si]", "3"},
  {"[Pb]", "3"},
  {"[Ge;!D5]", "3"},
  {"[Sn;!D5]", "3"},
  {"[AsD3]", "3"},
  {"[SbD3]", "3"},
  {"[BiD3]", "3"},
  {"[BiD5]", "6"},
  {"[se]", "2"},
  {"[Se]", "3"},
  {"[Te]", "3"},
  {"[Po]", "3"},
  {"[Be]", "1"},
  {"[Mg]", "1"},
  {"[Ca]", "1"},
  {"[Sr]", "1"},
  {"[Ba]", "1"},
  {"[Ra]", "1"}
};

class _Pattern {
public:

  _Pattern(string & format, int hyb_t) : smart_format(format), hyb_type(hyb_t) {
  }

  int get_type() {
    return this->hyb_type;
  }

  const string & get_smart_format() {
    return this->smart_format;
  }

  string smart_format;
  int hyb_type;
};

static void _parse(vector <_Pattern> & hybs) {
  for (const vector<string> &hyb_pattern : hyb_patterns) {
    if (hyb_pattern.size() != 2) {
      throw HessException("Incorrect representation of hybridization patterns");
    }
    string smart_buff = "";
    string hyb_buff = "";
    smart_buff = hyb_pattern[0];
    hyb_buff = hyb_pattern[1];
    int type = atoi(hyb_buff.c_str());
    if (type < 1 || type > 6) {
      throw HessException("Hyb types file invalid format. Incorrect hybridization");
    }
    _Pattern curr_hyb(smart_buff, type);
    hybs.push_back(curr_hyb);
  }
}

void hybridize(hess::Molecule* mol) {
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    mol->get_atom(a_id)->hybtype = 0;
  }
  vector <_Pattern> hyb_patterns;
  _parse(hyb_patterns);
  for (int index = 0; index < hyb_patterns.size(); index++) {
    _Pattern& pattern = hyb_patterns[index];
    const string & smarts_format = pattern.get_smart_format();
    const int& hyb_type = pattern.get_type();
    indigo::BufferScanner scanner(smarts_format.c_str());
    indigo::SmilesLoader loader(scanner);
    indigo::QueryMolecule qmol;
    loader.loadSMARTS(qmol);
    indigo::Molecule molcopy;
    indigo::Array<int> molcopy_mapping;
    molcopy.clone(*mol, &molcopy_mapping);
    indigo::AromaticityOptions arom_opt;
    molcopy.aromatize(arom_opt);
    indigo::MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qmol);
    indigo::MoleculeSubstructureMatcher::FragmentMatchCache fmcache;
    matcher.fmcache = &fmcache;

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[0]])->hybtype = hyb_type;
    }
  }
  // check all atoms to make sure *some* hybridization is assigned
  for (int a_id = mol->vertexBegin(); a_id != mol->vertexEnd(); a_id = mol->vertexNext(a_id)) {
    hess::Atom *a = mol->get_atom(a_id);
    if (a->hybtype == 0) {
      int ex_bond_count = mol->getVertex(a_id).degree();
      switch (ex_bond_count) {
        case 0:
        case 1:
        case 2:
          a->hybtype = 1;
          break;
        case 3:
          a->hybtype = 2;
          break;
        case 4:
          a->hybtype = 3;
          break;
        default:
          if (ex_bond_count == 5) {
            a->hybtype = 5;
          } else if (ex_bond_count == 6) {
            a->hybtype = 6;
          }
      }
    }
  }
}

void antialiasing(hess::Molecule* mol) {
  for (int atom_id = mol->vertexBegin(); atom_id != mol->vertexEnd(); atom_id = mol->vertexNext(atom_id)) {
    hess::Atom *atom = mol->get_atom(atom_id);
    if (mol->getVertex(atom_id).degree() < 2 || atom->hybtype == 0) {
      continue;
    }
    if (atom->hybtype == 1) {
      bool open_nb = false;
      ver atom_ver = mol->getVertex(atom_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        hess::Atom *b = mol->get_atom(nei_id);
        if ((b->hybtype == 1) || mol->getVertex(nei_id).degree() == 1) {
          open_nb = true;
        }
      }
      if (!open_nb) {
        atom->hybtype = 2;
      }
    } else if (atom->hybtype == 2) {
      bool open_nb = false;
      ver atom_ver = mol->getVertex(atom_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_id = atom_ver.neiVertex(id);
        hess::Atom *b = mol->get_atom(nei_id);
        if (b->hybtype == 2 || mol->getVertex(nei_id).degree() == 1) {
          open_nb = true;
        }
      }
      if (!open_nb) {
        atom->hybtype = 3;
      }
    }
  }//end 6c
}

double get_angle(const hess::Atom & a, const hess::Atom & b, const hess::Atom & c) {
  double ax, ay, az, cx, cy, cz, mdla, mdlb, mdlab, ccos, ugl;
  ax = b.x - a.x; //ab
  ay = b.y - a.y;
  az = b.z - a.z;
  cx = c.x - a.x; //ac
  cy = c.y - a.y;
  cz = c.z - a.z;
  int gr = 0, mn = 0, sc = 0;
  mdla = sqrt(ax * ax + ay * ay + az * az);
  mdlb = sqrt(cx * cx + cy * cy + cz * cz);
  mdlab = mdla * mdlb;
  ccos = (ax * cx + ay * cy + az * cz) / mdlab;
  gr = acos(ccos) * 180 / M_PI;
  return gr;
}

void hydridize_atoms(hess::Molecule* mol) {
  // step 6a
  // Set the levels of hybridization based on the average values ​​of the angles that the non-terminal atom makes with its neighbors.
  for (int id = mol->vertexBegin(); id != mol->vertexEnd(); id = mol->vertexNext(id)) {
    double sum_angle = 0;
    vector <int> neis_id;
    ver atom_ver = mol->getVertex(id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      neis_id.push_back(atom_ver.neiVertex(id));
    }
    if (neis_id.size() < 2) {
      continue;
    }
    set <pair<int, int>> exsp;
    hess::Atom *a = mol->get_atom(id);
    int bonds_count = 0;
    for (int b_id : neis_id) {
      for (int c_id : neis_id) {
        if (b_id == c_id) {
          continue;
        }
        if (exsp.find(make_pair(c_id, b_id)) != exsp.end()) {
          continue;
        }
        hess::Atom *b = mol->get_atom(b_id);
        hess::Atom *c = mol->get_atom(c_id);
        double ang = get_angle(*a, *b, *c);
        sum_angle += ang;
        bonds_count += 1;
        exsp.insert(make_pair(b_id, c_id));
      }
    }
    double av_angle = sum_angle / bonds_count;
    if (av_angle > 155) {
      a->hybtype = 1;
    } else if (av_angle > 115) {
      a->hybtype = 2;
    }      // (from obabel modify) -- special case for imines
    else if (mol->getAtomNumber(id) == 7 && neis_id.size() == 2 && av_angle > 109.5) {
      int h_count = 0;
      for (int b_id : neis_id) {
        if (mol->getAtomNumber(b_id) == 1) {
          h_count += 1;
        }
      }
      if (h_count == 1) {
        a->hybtype = 2;
      }
    } else {
      a->hybtype = 3;
    }
  }
}
