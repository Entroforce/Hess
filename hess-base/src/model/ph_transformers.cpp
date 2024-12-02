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

#include "bond/bond_order_determination.h"
#include "base_cpp/scanner.h"
#include "molecule/smiles_loader.h"
#include "molecule/query_molecule.h"
#include "reaction/rsmiles_loader.h"
#include "reaction/query_reaction.h"
#include "molecule/molecule_substructure_matcher.h"

using namespace hess;
using namespace indigo;

const vector<vector<string>> pH_data = {
  {"O=C[OD1-0:1]", "O=C[O-:1]", "4.0"},
  {"[O]=[C][C]=[C][O:1]", "[O]=[C][C]=[C][O-:1]", "4.0"},
  {"[N^3;!$(N~[!#6;!#1]):1]", "[N+:1]", "10.0"},
  {"[nD2:1]1c[nH]cc1", "[n+:1]1c[nH]cc1", "7.0"},
  {"[ND3+0:1]=[#6]", "[ND3+:1]=[#6]", "4.0"},
  {"[nD2:1]([#1])1[nD2-0][nD2-0][nD2-0]c1", "[n-:1]1nnnc1", "4.89"}, //?
  {"[nD2-0]1[nD2:1]([#1])[nD2-0][nD2-0]c1", "n1[n-:1]nnc1", "4.89"},
  {"[nD2-0:1]1[nD2-0][nD2-0][nD2-0]c1", "[n-:1]1nnnc1", "4.89"},
  {"[ND1:1]~[ND2:2]~[ND1:3]", "[N-:1]=[N+:2]=[N-:3]", "1E+10"},
  {"[ND1:1]~[ND2:2]~[ND1]~*", "[N-:1]=[N+:2]=[N]-*", "1E+10"},
  {"[NH2:1]~[NH1:2]~[NH1]~*", "[N-:1]=[N+:2]=[N]-*", "1E+10"},
  {"[NH1:1]~[ND2:2]~[NH1]~*", "[N-:1]=[N+:2]=[N]-*", "1E+10"},
  {"[O]=[N:1]-[O:2]", "[O]=[N+:1]-[O-:2]", "1E+10"},
  {"[O]-[N:1]-[O:2]", "[O]=[N+:1]-[O-:2]", "1E+10"},
  {"O=CN[OD1-0:1]", "O=CN[O-:1]", "8.0"},
  {"[SD3](=O)[OD1:1]", "[SD3](=O)[O-:1]", "2.0"},
  {"[SD3](=O)[O:1][#1]", "[SD3](=O)[O-:1]", "2.0"},
  {"[SD4]([!D1])(=O)(=O)[OD1:1]", "[SD4]([!D1])(=O)(=O)[O-:1]", "-2.6"},
  {"[SD4]([!D1])(=O)(=O)[O:1][#1]", "[SD4]([!D1])(=O)(=O)[O-:1]", "-2.6"},
  {"[#6^2+0](=[N^2+0:1])(~[N^2])*", "[#6](=[N+1:1])(~[N])*", "12.5"},
  {"[PD4](=O)([OD2])([OD2])[OD1:1]", "[PD4](=O)([OD2])([OD2])[O-:1]", "2.0"},
  {"[PD4](=O)([OD2])([!#8])[OD1:1]", "[PD4](=O)([OD2])([!#8])[O-:1]", "2.2"},
  //    {"O=P([!D1])([O:1][#1:2])[OD1]", "O=P([!D1])([O:1])O", "2.12"},
  {"O=P([*D2,*D3])([OD1:1])[OD1:2]", "O=P([*D2,*D3])([O-:1])[O-:2]", "2.12"},
  {"O=C(O)C(N)CC(=O)[OD1:1]", "O=C(O)C(N)CC(=O)[O-:1]", "3.8"},
  {"O=C(NCC=O)C(N)CC(=O)[OD1:1]", "O=C(NCC=O)C(N)CC(=O)[O-:1]", "3.8"},
  {"O=C(O)C(N)CCC(=O)[OD1:1]", "O=C(O)C(N)CCC(=O)[O-:1]", "5.0"},
  {"O=C(NCC=O)C(N)CCC(=O)[OD1:1]", "O=C(NCC=O)C(N)CCC(=O)[O-:1]", "5.0"},
  {"O=C(O)C(N)CCCNC(N)=[N:1]", "O=C(O)C(N)CCCNC(N)=[N+:1]", "12.0"},
  {"O=C(NCC=O)C(N)CCCNC(N)=[N:1]", "O=C(NCC=O)C(N)CCCNC(N)=[N+:1]", "12.0"},
  {"O=C(O)C(N)CCCC[N:1]", "O=C(O)C(N)CCCC[N+:1]", "10.5"},
  {"O=C(NCC=O)C(N)CCCC[N:1]", "O=C(NCC=O)C(N)CCCC[N+:1]", "10.5"},
  {"O=C(O)C(N)Cc1[nH0:1]c[nH]c1", "O=C(O)C(N)Cc1[n+:1]c[nH]c1", "6.08"},
  {"O=C(O)C(N)Cc1[nH]c[nH0:1]c1", "O=C(O)C(N)Cc1[nH]c[n+:1]c1", "6.08"},
  {"O=C(NCC=O)C(N)Cc1[nH0:1]c[nH]c1", "O=C(NCC=O)C(N)Cc1[n+:1]c[nH]c1", "6.08"},
  {"O=C(O)C(N)C[S:1]", "O=C(O)C(N)C[S-:1]", "8.28"},
  {"O=C(NCC=O)C(N)C[S:1]", "O=C(NCC=O)C(N)C[S-:1]", "8.28"},
  {"O=C(O)C(N)Cc1ccc([O:1])cc1", "O=C(O)C(N)Cc1ccc([O-:1])cc1", "10.1"},
  {"O=C(NCC=O)C(N)Cc1ccc([O:1])cc1", "O=C(NCC=O)C(N)Cc1ccc([O-:1])cc1", "10.1"},
};

class PhTransform {
public:
  PhTransform(const std::string &reag, const std::string &prod, const double pKa);
  std::string reagents;
  std::string products;
  double pKa;
};

PhTransform::PhTransform(const string &reag, const string &prod, const double pKa) {
  this->pKa = pKa;
  this->reagents = reag;
  this->products = prod;
}

static void _parse(vector<PhTransform> & transformations) {
  for (const vector<string> &ph_tr : pH_data) {
    if (ph_tr.size() != 3) {
      throw HessException("Incorrect representation of pH transformation.\n");
    }
    double pKa = atof(ph_tr[2].c_str());
    PhTransform ptrans(ph_tr[0], ph_tr[1], pKa);
    transformations.push_back(ptrans);
  }
}

void transform_Ph(hess::Molecule* mol) {
  vector<PhTransform> transformations;
  _parse(transformations);

  indigo::Molecule molcopy;
  Array<int> molcopy_mapping;
  molcopy.clone(*mol, &molcopy_mapping);
  indigo::AromaticityOptions arom_options(indigo::AromaticityOptions::GENERIC);
  molcopy.aromatize(arom_options);

  for (PhTransform &tr : transformations) {

    std::vector<int> hyb_indices;
    std::vector<int> hyb_values;
    if (tr.reagents == "[N^3;!$(N~[!#6;!#1]):1]") {
      hyb_indices.push_back(0);
      hyb_values.push_back(3);
      tr.reagents = "[N;!$(N~[!#6;!#1]):1]";
    } else if (tr.reagents == "[#6^2+0](=[N^2+0:1])(~[N^2])*") {
      hyb_indices.push_back(0);
      hyb_indices.push_back(1);
      hyb_indices.push_back(2);
      hyb_values.push_back(2);
      hyb_values.push_back(2);
      hyb_values.push_back(2);
      tr.reagents = "[#6+0](=[N+0:1])(~[N])*";
    }

    indigo::QueryMolecule qreactant, qproduct;
    {
      indigo::BufferScanner scanner(tr.reagents.c_str());
      indigo::SmilesLoader loader(scanner);
      loader.inside_rsmiles = true; // to get AAM
      loader.loadSMARTS(qreactant);
    }
    {
      indigo::BufferScanner scanner(tr.products.c_str());
      indigo::SmilesLoader loader(scanner);
      loader.inside_rsmiles = true; // to get AAM
      loader.loadSMARTS(qproduct);
    }

    int totalchg = 0;
    for (int j = qproduct.vertexBegin(); j != qproduct.vertexEnd(); j = qproduct.vertexNext(j)) {
      int charge = qproduct.getAtomCharge(j);
      if (charge != indigo::CHARGE_UNKNOWN)
        totalchg += charge;
    }

    MoleculeSubstructureMatcher matcher(molcopy);
    matcher.setQuery(qreactant);
    MoleculeSubstructureMatcher::FragmentMatchCache fmcache;
    matcher.fmcache = &fmcache;

    for (bool flag = matcher.find(); flag; flag = matcher.findNext()) {
      int i;

      for (i = 0; i < hyb_indices.size(); i++)
        if (mol->get_atom(molcopy_mapping[matcher.getQueryMapping()[hyb_indices[i]]])->hybtype != hyb_values[i])
          break;
      if (i != hyb_indices.size())
        continue;
      
      if (tr.pKa < 1e9) {
        if (totalchg > 0 && pow(10, tr.pKa - mol->pH) < 1.0)
          continue;
        if (totalchg < 0 && pow(10, tr.pKa - mol->pH) > 1.0)
          continue;
      }

      for (i = qreactant.vertexBegin(); i != qreactant.vertexEnd(); i = qreactant.vertexNext(i)) {
        int aam = qreactant.getAAMArray()[i];
        if (aam > 0)
          for (int j = qproduct.vertexBegin(); j != qproduct.vertexEnd(); j = qproduct.vertexNext(j))
            if (qproduct.getAAMArray()[j] == aam) {
              int chg = qproduct.getAtomCharge(j);
              if (chg != indigo::CHARGE_UNKNOWN){
                mol->setAtomCharge(molcopy_mapping[matcher.getQueryMapping()[i]], chg);
              }
            }
      }
      
      if (qproduct.edgeCount() != qreactant.edgeCount()) // case [nD2:1]([#1])1[nD2-0][nD2-0][nD2-0]c1 >> [n-:1]1nnnc1 for CHEMBL lingand
        continue;
      
      for (i = qreactant.edgeBegin(); i != qreactant.edgeEnd(); i = qreactant.edgeNext(i)) {
        int order = qproduct.getBondOrder(i);
        if (order > 0){
           bond_order_changed(
                   molcopy_mapping[matcher.getQueryMapping()[qreactant.getEdge(i).beg]],
                   molcopy_mapping[matcher.getQueryMapping()[qreactant.getEdge(i).end]],
                   order, mol);
        }
      }
    }
  }
}
