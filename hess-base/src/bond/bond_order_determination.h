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

#include "model/molecule.h"
#include "model/vector3.h"

//! \return Returns the bond length in angstroms.
double distance(const hess::Atom & a, const hess::Atom & b);

double distance_sqr(const hess::Atom & a, const hess::Atom & b);

void determine_connectivity(hess::Molecule* mol);

void delete_bonds_exceeding_valence(hess::Molecule* mol);

void bond_order_changed(const int a_id, const int b_id, const int bond_order, hess::Molecule* mol);

bool has_bond_of_order(int a_id, int bond_order, hess::Molecule* mol);

//! Assign double or triple bond based on hybridization and geometric features .
void assign_bond_orders(hess::Molecule* mol, bool after_aromatize = true);

void unfold_hydrogens(hess::Molecule* mol);

int count_explicit_h(int a_id, hess::Molecule* mol);

bool has_nonsp2(vector <int> & ring_atoms, hess::Molecule* mol);

int nonhydrobond_count(hess::Molecule* mol, int a_id);

//! \return The number of oxygen atoms connected that only have one heavy valence
int count_free_oxygens(hess::Molecule* mol, const int a_id);

//! \return Is this atom an oxygen in a carboxyl (-CO2 or CO2H) group
bool is_carboxyl_oxygen(hess::Molecule* mol, int a_id);


//! Assign hybridization for atoms based on the values ​​of the average angles 
void hydridize_atoms(hess::Molecule* mol);
//! \return angle between the vectors formed by the coordinates of atoms a,b and a,c
double get_angle(const hess::Atom & a, const hess::Atom & b, const hess::Atom & c);
//! Assign new orders of hybridization - reduce them depending on the neighbors of the atom.
void antialiasing(hess::Molecule* mol);

//! Assign all ring compounds in a molecule
void extract_rings(vector <vector <int>> &rings_bonds, vector <vector <int>> &rings_atoms, hess::Molecule *mol);
//! Assign sp2 hybridization for 5-6 ring atoms depending on the value of the twist angle
void check_torsion_angles(hess::Molecule* mol, vector <vector <int>> &rings_atoms);
//! \return twist angle for atoms a b c d
double get_torsion_angle(int a_num, int b_num, int c_num, int d_num, hess::Molecule* mol);

void determine_bond_orders_in_aromatic_rings(hess::Molecule* hess_mol, vector <vector <int>> &rings_bonds, vector <vector <int>> &rings_atoms);


//! Assign functional groups using SMARTS templates
void functional_groups(hess::Molecule* mol);

void hybridize(hess::Molecule* mol);

void transform_Ph(hess::Molecule* mol);
