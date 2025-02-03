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

/* Some parts of this file were borrowed from
 * https://github.com/openbabel/openbabel/blob/master/src/formats/pdbformat.cpp
 * Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
 * Some portions Copyright (C) 2003-2006 Geoffrey R. Hutchison
 * Some portions Copyright (C) 2004 by Chris Morley
 * Part of the Open Babel package, under the GNU General Public License (GPL)
 */


#include "parser/parser.h"
#include "molecule/molecule_auto_loader.h"
#include "base_cpp/scanner.h"
#include "elements.h"

#include <fstream>
using namespace hess;
using namespace std;

/*

  From http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

     COLUMNS        DATA TYPE       CONTENTS
     --------------------------------------------------------------------------------
      1 -  6        Record name     "ATOM  "
      7 - 11        Integer         Atom serial number.
     13 - 16        Atom            Atom name.
     17             Character       Alternate location indicator.
     18 - 20        Residue name    Residue name.
     22             Character       Chain identifier.
     23 - 26        Integer         Residue sequence number.
     27             AChar           Code for insertion of residues.
     31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
     39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
     47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
     55 - 60        Real(6.2)       Occupancy.
     61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
     73 - 76        LString(4)      Segment identifier, left-justified.
     77 - 78        LString(2)      Element symbol, right-justified.
     79 - 80        LString(2)      Charge on the atom.
 */



void del_tabs(string& in_str) {
  for (int i = 0; i < in_str.size(); i++) {
    if (in_str[i] == ' ' || in_str[i] == '\t' || in_str[i] == '\r') {
      in_str.erase(i, 1);
      i--;
    }
  }
}

unsigned int parse_element_from_type(string &sbuf, bool hetatm) {
  /* atom name */
  string atmid = sbuf.substr(6, 4);

  /* element */
  string element = "  ";


  //trim spaces on the right and left sides
  while (!atmid.empty() && atmid[0] == ' ')
    atmid = atmid.erase(0, 1);

  while (!atmid.empty() && atmid[atmid.size() - 1] == ' ')
    atmid = atmid.substr(0, atmid.size() - 1);

  /* residue name */
  string resname = sbuf.substr(11, 3);
  if (resname == "   ")
    resname = "UNK";
  else {
    while (!resname.empty() && resname[0] == ' ')
      resname = resname.substr(1, resname.size() - 1);

    while (!resname.empty() && resname[resname.size() - 1] == ' ')
      resname = resname.substr(0, resname.size() - 1);
  }

  string type;
  // OK, we have to fall back to determining the element from the atom type
  // This is unreliable, but there's no other choice
  if (!hetatm) {
    type = atmid.substr(0, 2);
    if (isdigit(type[0])) {
      // sometimes non-standard files have, e.g 11HH
      if (!isdigit(type[1])) type = atmid.substr(1, 1);
      else type = atmid.substr(2, 1);
    }
    else if ((sbuf[6] == ' ' &&
            strncasecmp(type.c_str(), "Zn", 2) != 0 &&
            strncasecmp(type.c_str(), "Fe", 2) != 0) ||
            isdigit(type[1])) //type[1] is digit in Platon
      type = atmid.substr(0, 1); // one-character element


    if (resname.substr(0, 2) == "AS" || resname[0] == 'N') {
      if (atmid == "AD1")
        type = "O";
      if (atmid == "AD2")
        type = "N";
    }
    if (resname.substr(0, 3) == "HIS" || resname[0] == 'H') {
      if (atmid == "AD1" || atmid == "AE2")
        type = "N";
      if (atmid == "AE1" || atmid == "AD2")
        type = "C";
    }
    if (resname.substr(0, 2) == "GL" || resname[0] == 'Q') {
      if (atmid == "AE1")
        type = "O";
      if (atmid == "AE2")
        type = "N";
    }
    // fix: #2002557
    if (atmid[0] == 'H' &&
            (atmid[1] == 'D' || atmid[1] == 'E' ||
            atmid[1] == 'G' || atmid[1] == 'H' ||
            atmid[1] == 'N')) // HD, HE, HG, HH, HN...
      type = "H";

    if (type.size() == 2)
      type[1] = tolower(type[1]);

  }
  else { //must be hetatm record
    if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' '))) {
      if (isalpha(element[0]))
        type = element.substr(0, 2);
      else
        type = element.substr(1, 1);

      if (type.size() == 2)
        type[1] = tolower(type[1]);
    }
    else { // no element column to use
      if (isalpha(atmid[0])) {
        if (atmid.size() > 2)
          type = atmid.substr(0, 2);
        else if (atmid[0] == 'A') // alpha prefix
          type = atmid.substr(1, atmid.size() - 1);
        else
          type = atmid.substr(0, 1);
      }
      else if (atmid[0] == ' ')
        type = atmid.substr(1, 1); // one char element
      else
        type = atmid.substr(1, 2);

      // Some cleanup steps
      if (atmid == resname) {
        type = atmid;
        if (type.size() == 2)
          type[1] = tolower(type[1]);
      }
      else
        if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
              resname == "GPG" || resname == "NAD" || resname == "NAL" ||
              resname == "NDP" || resname == "ABA") {
        if (type.size() > 1)
          type = type.substr(0, 1);
        //type.erase(1,type.size()-1);
      }
      else // other residues
        if (isdigit(type[0])) {
        type = type.substr(1, 1);
      } else if (type.size() > 1 && isdigit(type[1]))
        type = type.substr(0, 1);
      else if (type.size() > 1 && isalpha(type[1])) {
        if (type[0] == 'O' && type[1] == 'H')
          type = type.substr(0, 1); // no "Oh" element (e.g. 1MBN)
        else if (isupper(type[1])) {
          type[1] = tolower(type[1]);
        }
      }
    }

  } // HETATM records

  unsigned int atomic_num = GetAtomicNum(type.c_str());
  if (atomic_num == 0) { //try one character if two character element not found
    type = type.substr(0, 1);
    atomic_num = GetAtomicNum(type.c_str());
  }
  //    atom.SetAtomicNum(atomic_num);
  return atomic_num;

}

bool is_mol_format(const char* file_name, ifstream& ifile) {
  const string& name_buf = string(file_name);
  int size = name_buf.size();
  if (name_buf.substr(size - 3, 3) != "mol")
    return false;
  string curr_str;
  int str_num = 1;

  while (getline(ifile, curr_str)) {
    if (str_num == 4) {
      del_tabs(curr_str);
      const string& buf = curr_str.substr(curr_str.size() - 5, 5);
      if (buf != "V3000" && buf != "V2000")
        return false;
      break;
    }
    str_num++;
  }
  return true;
}

void parse_pdb(hess::Molecule* molecule, const char* file_name, ifstream& file) {
  string curr_str;
  int row_i = 1;
  int check = 0;
  int offset = 0;
  while (getline(file, curr_str)) {
    string type = curr_str.substr(0, 6);
    del_tabs(type);
    if (type == "ATOM" || type == "HETATM") {
        if(check == 0){
            check = 1;
            offset = row_i - 1;
        }
      char *buffer = const_cast<char*> (curr_str.c_str());
      string sbuf = &buffer[6];
      bool hetatm = (type == "HETATM") ? true : false;
      bool elementFound = false; // true if correct element found in col 77-78
      if (sbuf.size() < 48) {
        throw HessException("In " + string(file_name) + " invalid pdb file format. Check line number " + to_string(row_i));
      }
      string atom_name = curr_str.substr(12, 4);
      del_tabs(atom_name);
      string residue_name = curr_str.substr(17, 3);
      del_tabs(residue_name);
      string res_seq_number = curr_str.substr(22, 4);
      del_tabs(res_seq_number);
      int res_seq = stoi(res_seq_number);
      string x_str = curr_str.substr(30, 8);
      del_tabs(x_str);
      double x = stod(x_str);
      string y_str = curr_str.substr(38, 8);
      del_tabs(y_str);
      double y = stod(y_str);
      string z_str = curr_str.substr(46, 8);
      del_tabs(z_str);
      double z = stod(z_str);
      string element = "  ";
      unsigned int element_num = -1;
      if (sbuf.size() > 71) {
        element = sbuf.substr(70, 2);
        if (isalpha(element[1])) {
          if (element[0] == ' ') {
            element.erase(0, 1);
            elementFound = true;
          } else if (isalpha(element[0])) {
            elementFound = true;
            element[1] = tolower(element[1]);
          }
        }
      }
      if (!elementFound) {
        element_num = parse_element_from_type(sbuf, hetatm);
      } else {
        element_num = GetAtomicNum(element.c_str());
      }

      if (element_num == -1 || element_num == 0) {
        if (atom_name.size() > 1)
          atom_name[1] = std::tolower(atom_name[1]);
        element_num = GetAtomicNum(atom_name.c_str());
        if (element_num == -1 || element_num == 0)
          throw HessException("Failed to recognize element. Check line number " + to_string(row_i));
      }
      int charge = 0;
      if (sbuf.size() > 73) {
        string charge_str = curr_str.substr(78, 2);
        del_tabs(charge_str);
//        std::cout << row_i << " " << charge_str << std::endl;
        if (charge_str != "") {
          charge = charge_str[0] - '0';
          if (charge_str[1] == '-') {
            charge *= -1;
          }
        }
      }
      
      char txtatmi[26] = {0};
//      snprintf(txtatmi, 26, "( %5d %3s %4d)        ", row_i - offset, residue_name.c_str(), res_seq);
//      std::cout << txtatmi << " --- " << residue_name << " " << res_seq << std::endl;
      ((Molecule*) molecule)->txtatm = (char **) realloc(((Molecule*) molecule)->txtatm, (row_i - offset) * sizeof(char *));
      ((Molecule*) molecule)->txtatm[row_i - offset - 1] = (char *) malloc(26 * sizeof(char));
      
      memcpy(((Molecule*) molecule)->txtatm[row_i - offset - 1], txtatmi, 26);
      
      int id = ((Molecule*) molecule)->add_atom(x, y, z, element_num);
      ((Molecule*) molecule)->setAtomCharge(id, charge);
    }
    row_i++;
  }
}

void parse_mol(hess::Molecule* molecule, const char* file_name) {
  hess::Molecule &mol = *(hess::Molecule *)molecule;
  indigo::FileScanner scanner(file_name);
  indigo::MoleculeAutoLoader loader(scanner);

  loader.loadMolecule(mol);

  mol.atoms.resize(mol.vertexEnd());
  mol.coords.resize(mol.vertexEnd());

  for (int i = mol.vertexBegin(); i != mol.vertexEnd(); i = mol.vertexNext(i)) {
    mol.atoms[i].x = mol.getAtomXyz(i).x;
    mol.atoms[i].y = mol.getAtomXyz(i).y;
    mol.atoms[i].z = mol.getAtomXyz(i).z;

    mol.atoms[i].x_orign = mol.getAtomXyz(i).x;
    mol.atoms[i].y_orign = mol.getAtomXyz(i).y;
    mol.atoms[i].z_orign = mol.getAtomXyz(i).z;

    mol.coords[i].x = mol.getAtomXyz(i).x;
    mol.coords[i].y = mol.getAtomXyz(i).y;
    mol.coords[i].z = mol.getAtomXyz(i).z;
  }
}

hess::Molecule* Parser::parse_molecule(const char* file_name) {
  if (!file_name) {
    throw HessException("The path to the molecule is not set");
  }
  ifstream file;
  hess::Molecule* molecule = new Molecule();
  file.open(file_name);
  if (!file) {
    throw HessException(string(file_name) + " file open error");
  }
  try {
    if (is_mol_format(file_name, file)) {
      parse_mol(molecule, file_name);
    } else {
      parse_pdb(molecule, file_name, file);
    }
    if (!((Molecule*) molecule)->get_atoms_count()) {
      throw HessException("No atoms in " + string(file_name));
    }
  } catch (HessException &e) {
    delete molecule;
    throw HessException(e.what());
  } catch (indigo::Exception &e) {
    delete molecule;
    throw HessException(e.what());
  } catch (exception &e) {
    delete molecule;
    throw HessException(e.what());
  }
  return molecule;
}
