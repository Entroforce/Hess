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
 * https://github.com/mwojcikowski/smina/blob/master/src/lib/atom_constants.h
 * published under Apache License 2.0
 * Copyright (c) 2006-2010, The Scripps Research Institute 
 ***********************************************************************/

#pragma once

#include <string>
#include <assert.h>

// based on SY_TYPE_* but includes H
const size_t EL_TYPE_H = 0;
const size_t EL_TYPE_C = 1;
const size_t EL_TYPE_N = 2;
const size_t EL_TYPE_O = 3;
const size_t EL_TYPE_S = 4;
const size_t EL_TYPE_P = 5;
const size_t EL_TYPE_F = 6;
const size_t EL_TYPE_Cl = 7;
const size_t EL_TYPE_Br = 8;
const size_t EL_TYPE_I = 9;
const size_t EL_TYPE_Met = 10;
const size_t EL_TYPE_SIZE = 11;

// AutoDock4
const size_t AD_TYPE_C = 0;
const size_t AD_TYPE_A = 1;
const size_t AD_TYPE_N = 2;
const size_t AD_TYPE_O = 3;
const size_t AD_TYPE_P = 4;
const size_t AD_TYPE_S = 5;
const size_t AD_TYPE_H = 6; // non-polar hydrogen
const size_t AD_TYPE_F = 7;
const size_t AD_TYPE_I = 8;
const size_t AD_TYPE_NA = 9;
const size_t AD_TYPE_OA = 10;
const size_t AD_TYPE_SA = 11;
const size_t AD_TYPE_HD = 12;
const size_t AD_TYPE_Mg = 13;
const size_t AD_TYPE_Mn = 14;
const size_t AD_TYPE_Zn = 15;
const size_t AD_TYPE_Ca = 16;
const size_t AD_TYPE_Fe = 17;
const size_t AD_TYPE_Cl = 18;
const size_t AD_TYPE_Br = 19;
const size_t AD_TYPE_METAL = 20; //generic metal, not actually part of autodock
const size_t AD_TYPE_SIZE = 21;

// X-Score
const size_t XS_TYPE_C_H = 0;
const size_t XS_TYPE_C_P = 1;
const size_t XS_TYPE_N_P = 2;
const size_t XS_TYPE_N_D = 3;
const size_t XS_TYPE_N_A = 4;
const size_t XS_TYPE_N_DA = 5;
const size_t XS_TYPE_O_P = 6;
const size_t XS_TYPE_O_D = 7;
const size_t XS_TYPE_O_A = 8;
const size_t XS_TYPE_O_DA = 9;
const size_t XS_TYPE_S_P = 10;
const size_t XS_TYPE_P_P = 11;
const size_t XS_TYPE_F_H = 12;
const size_t XS_TYPE_Cl_H = 13;
const size_t XS_TYPE_Br_H = 14;
const size_t XS_TYPE_I_H = 15;
const size_t XS_TYPE_Met_D = 16;
const size_t XS_TYPE_SIZE = 17;


namespace specific_atom_type {

  enum type {
    Hydrogen, // H_H_X,
    PolarHydrogen, //(can donate) H_HD_X,
    AliphaticCarbonXSHydrophobe, // C_C_C_H, //hydrophobic according to xscale
    AliphaticCarbonXSNonHydrophobe, //C_C_C_P, //not hydrophobic (according to xs)
    AromaticCarbonXSHydrophobe, //C_A_C_H,
    AromaticCarbonXSNonHydrophobe, //C_A_C_P,
    Nitrogen, //N_N_N_P, no hydrogen bonding
    NitrogenXSDonor, //N_N_N_D,
    NitrogenXSDonorAcceptor, //N_NA_N_DA, also an autodock acceptor
    NitrogenXSAcceptor, //N_NA_N_A, also considered an acceptor by autodock
    Oxygen, //O_O_O_P,
    OxygenXSDonor, //O_O_O_D,
    OxygenXSDonorAcceptor, //O_OA_O_DA, also an autodock acceptor
    OxygenXSAcceptor, //O_OA_O_A, also an autodock acceptor
    Sulfur, //S_S_S_P,
    SulfurAcceptor, //S_SA_S_P, XS doesn't do sulfur acceptors
    Phosphorus, //P_P_P_P,
    Fluorine, //F_F_F_H,
    Chlorine, //Cl_Cl_Cl_H,
    Bromine, //Br_Br_Br_H,
    Iodine, //I_I_I_H,
    Magnesium, //Met_Mg_Met_D,
    Manganese, //Met_Mn_Met_D,
    Zinc, // Met_Zn_Met_D,
    Calcium, //Met_Ca_Met_D,
    Iron, //Met_Fe_Met_D,
    GenericMetal, //Met_METAL_Met_D,
    NumTypes
  };

  //store all the desired properties in atom_type_info

  struct info {
    type sm;
    size_t el;
    size_t ad;
    size_t xs;
    const char* special_name; //this must be more than 2 chars long
    const char* adname; //this must be no longer than 2 chars
    double ad_radius;
    double ad_depth;
    double ad_solvation;
    double ad_volume;
    double covalent_radius;
    double xs_radius;
    bool xs_hydrophobe;
    bool xs_donor;
    bool xs_acceptor;
    bool ad_heteroatom;

  };

  extern info data[NumTypes];


  const info default_data[NumTypes] = {//el, ad, xs
    {Hydrogen, EL_TYPE_H, AD_TYPE_H, XS_TYPE_SIZE, "Hydrogen",
      "H", 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false, false, false, false},
    {PolarHydrogen, EL_TYPE_H, AD_TYPE_HD, XS_TYPE_SIZE, "PolarHydrogen",
      "HD", 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false, false, false, false},
    {AliphaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, "AliphaticCarbonXSHydrophobe",
      "C", 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 1.900000, true, false, false, false},
    {AliphaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_P, "AliphaticCarbonXSNonHydrophobe",
      "C", 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 1.900000, false, false, false, false},
    {AromaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_H, "AromaticCarbonXSHydrophobe",
      "A", 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, true, false, false, false},
    {AromaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_P, "AromaticCarbonXSNonHydrophobe",
      "A", 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, false, false, false, false},
    {Nitrogen, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_P, "Nitrogen",
      "N", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, false, false, true},
    {NitrogenXSDonor, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_D, "NitrogenXSDonor",
      "N", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, true, false, true},
    {NitrogenXSDonorAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_DA, "NitrogenXSDonorAcceptor",
      "NA", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, true, true, true},
    {NitrogenXSAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_A, "NitrogenXSAcceptor",
      "NA", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, false, true, true},
    {Oxygen, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_P, "Oxygen",
      "O", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, false, false, true},
    {OxygenXSDonor, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_D, "OxygenXSDonor",
      "O", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, true, false, true},
    {OxygenXSDonorAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_DA, "OxygenXSDonorAcceptor",
      "OA", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, true, true, true},
    {OxygenXSAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_A, "OxygenXSAcceptor",
      "OA", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, false, true, true},
    {Sulfur, EL_TYPE_S, AD_TYPE_S, XS_TYPE_S_P, "Sulfur",
      "S", 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, false, false, false, true},
    {SulfurAcceptor, EL_TYPE_S, AD_TYPE_SA, XS_TYPE_S_P, "SulfurAcceptor",
      "SA", 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, false, false, false, true},
    {Phosphorus, EL_TYPE_P, AD_TYPE_P, XS_TYPE_P_P, "Phosphorus",
      "P", 2.100000, 0.200000, -0.001100, 38.792400, 1.060000, 2.100000, false, false, false, true},
    {Fluorine, EL_TYPE_F, AD_TYPE_F, XS_TYPE_F_H, "Fluorine",
      "F", 1.545000, 0.080000, -0.001100, 15.448000, 0.710000, 1.500000, true, false, false, true},
    {Chlorine, EL_TYPE_Cl, AD_TYPE_Cl, XS_TYPE_Cl_H, "Chlorine",
      "Cl", 2.045000, 0.276000, -0.001100, 35.823500, 0.990000, 1.800000, true, false, false, true},
    {Bromine, EL_TYPE_Br, AD_TYPE_Br, XS_TYPE_Br_H, "Bromine",
      "Br", 2.165000, 0.389000, -0.001100, 42.566100, 1.140000, 2.000000, true, false, false, true},
    {Iodine, EL_TYPE_I, AD_TYPE_I, XS_TYPE_I_H, "Iodine",
      "I", 2.360000, 0.550000, -0.001100, 55.058500, 1.330000, 2.200000, true, false, false, true},
    {Magnesium, EL_TYPE_Met, AD_TYPE_Mg, XS_TYPE_Met_D, "Magnesium",
      "Mg", 0.650000, 0.875000, -0.001100, 1.560000, 1.300000, 1.200000, false, true, false, true},
    {Manganese, EL_TYPE_Met, AD_TYPE_Mn, XS_TYPE_Met_D, "Manganese",
      "Mn", 0.650000, 0.875000, -0.001100, 2.140000, 1.390000, 1.200000, false, true, false, true},
    {Zinc, EL_TYPE_Met, AD_TYPE_Zn, XS_TYPE_Met_D, "Zinc",
      "Zn", 0.740000, 0.550000, -0.001100, 1.700000, 1.310000, 1.200000, false, true, false, true},
    {Calcium, EL_TYPE_Met, AD_TYPE_Ca, XS_TYPE_Met_D, "Calcium",
      "Ca", 0.990000, 0.550000, -0.001100, 2.770000, 1.740000, 1.200000, false, true, false, true},
    {Iron, EL_TYPE_Met, AD_TYPE_Fe, XS_TYPE_Met_D, "Iron",
      "Fe", 0.650000, 0.010000, -0.001100, 1.840000, 1.250000, 1.200000, false, true, false, true},
    {GenericMetal, EL_TYPE_Met, AD_TYPE_METAL, XS_TYPE_Met_D, "GenericMetal",
      "M", 1.200000, 0.000000, -0.001100, 22.449300, 1.750000, 1.200000, false, true, false, true}
  };

  const info vinardo_data[NumTypes] = {//el, ad, xs
    {Hydrogen, EL_TYPE_H, AD_TYPE_H, XS_TYPE_SIZE, "Hydrogen",
      "H", 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false, false, false, false},
    {PolarHydrogen, EL_TYPE_H, AD_TYPE_HD, XS_TYPE_SIZE, "PolarHydrogen",
      "HD", 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false, false, false, false},
    {AliphaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H, "AliphaticCarbonXSHydrophobe",
      "C", 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 2.000000, true, false, false, false},
    {AliphaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_P, "AliphaticCarbonXSNonHydrophobe",
      "C", 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 2.000000, false, false, false, false},
    {AromaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_H, "AromaticCarbonXSHydrophobe",
      "A", 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, true, false, false, false},
    {AromaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_P, "AromaticCarbonXSNonHydrophobe",
      "A", 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, true, false, false, false},
    {Nitrogen, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_P, "Nitrogen",
      "N", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.700000, false, false, false, true},
    {NitrogenXSDonor, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_D, "NitrogenXSDonor",
      "N", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.700000, false, true, false, true},
    {NitrogenXSDonorAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_DA, "NitrogenXSDonorAcceptor",
      "NA", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.700000, false, true, true, true},
    {NitrogenXSAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_A, "NitrogenXSAcceptor",
      "NA", 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.700000, false, false, true, true},
    {Oxygen, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_P, "Oxygen",
      "O", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.600000, false, false, false, true},
    {OxygenXSDonor, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_D, "OxygenXSDonor",
      "O", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.600000, false, true, false, true},
    {OxygenXSDonorAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_DA, "OxygenXSDonorAcceptor",
      "OA", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.600000, false, true, true, true},
    {OxygenXSAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_A, "OxygenXSAcceptor",
      "OA", 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.600000, false, false, true, true},
    {Sulfur, EL_TYPE_S, AD_TYPE_S, XS_TYPE_S_P, "Sulfur",
      "S", 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, false, false, false, true},
    {SulfurAcceptor, EL_TYPE_S, AD_TYPE_SA, XS_TYPE_S_P, "SulfurAcceptor",
      "SA", 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, true, false, false, true},
    {Phosphorus, EL_TYPE_P, AD_TYPE_P, XS_TYPE_P_P, "Phosphorus",
      "P", 2.100000, 0.200000, -0.001100, 38.792400, 1.060000, 2.100000, false, false, false, true},
    {Fluorine, EL_TYPE_F, AD_TYPE_F, XS_TYPE_F_H, "Fluorine",
      "F", 1.545000, 0.080000, -0.001100, 15.448000, 0.710000, 1.500000, true, false, false, true},
    {Chlorine, EL_TYPE_Cl, AD_TYPE_Cl, XS_TYPE_Cl_H, "Chlorine",
      "Cl", 2.045000, 0.276000, -0.001100, 35.823500, 0.990000, 1.800000, true, false, false, true},
    {Bromine, EL_TYPE_Br, AD_TYPE_Br, XS_TYPE_Br_H, "Bromine",
      "Br", 2.165000, 0.389000, -0.001100, 42.566100, 1.140000, 2.000000, true, false, false, true},
    {Iodine, EL_TYPE_I, AD_TYPE_I, XS_TYPE_I_H, "Iodine",
      "I", 2.360000, 0.550000, -0.001100, 55.058500, 1.330000, 2.200000, true, false, false, true},
    {Magnesium, EL_TYPE_Met, AD_TYPE_Mg, XS_TYPE_Met_D, "Magnesium",
      "Mg", 0.650000, 0.875000, -0.001100, 1.560000, 1.300000, 1.200000, false, true, false, true},
    {Manganese, EL_TYPE_Met, AD_TYPE_Mn, XS_TYPE_Met_D, "Manganese",
      "Mn", 0.650000, 0.875000, -0.001100, 2.140000, 1.390000, 1.200000, false, true, false, true},
    {Zinc, EL_TYPE_Met, AD_TYPE_Zn, XS_TYPE_Met_D, "Zinc",
      "Zn", 0.740000, 0.550000, -0.001100, 1.700000, 1.310000, 1.200000, false, true, false, true},
    {Calcium, EL_TYPE_Met, AD_TYPE_Ca, XS_TYPE_Met_D, "Calcium",
      "Ca", 0.990000, 0.550000, -0.001100, 2.770000, 1.740000, 1.200000, false, true, false, true},
    {Iron, EL_TYPE_Met, AD_TYPE_Fe, XS_TYPE_Met_D, "Iron",
      "Fe", 0.650000, 0.010000, -0.001100, 1.840000, 1.250000, 1.200000, false, true, false, true},
    {GenericMetal, EL_TYPE_Met, AD_TYPE_METAL, XS_TYPE_Met_D, "GenericMetal",
      "M", 1.200000, 0.000000, -0.001100, 22.449300, 1.750000, 1.200000, false, true, false, true}
  };
}

typedef specific_atom_type::type smt;

struct atom_equivalence {
  std::string name;
  std::string to;
};


const atom_equivalence atom_equivalence_data[] = {
  {"Se", "S"}
};


const size_t atom_equivalences_size = sizeof (atom_equivalence_data) / sizeof (const atom_equivalence);

const double epsilon_fl = std::numeric_limits<double>::epsilon();
const double max_fl = (std::numeric_limits<double>::max)();

const std::string non_ad_metal_names[] = {// expand as necessary
  "Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

inline bool is_non_ad_metal_name(const std::string& name) {
  const size_t s = sizeof (non_ad_metal_names) / sizeof (const std::string);
  for (int i = 0; i < s; i++) {
    if (non_ad_metal_names[i] == name)
      return true;
  }
  return false;
}

inline smt string_to_atom_type(const std::string& name) {
  if (name.length() <= 2) {
    for (int i = 0; i < smt::NumTypes; i++) {
      if (specific_atom_type::data[i].adname == name) {
        return specific_atom_type::data[i].sm;
      }
    }
    for (int i = 0; i < atom_equivalences_size; i++) {
      if (atom_equivalence_data[i].name == name)
        return string_to_atom_type(atom_equivalence_data[i].to);
    }
    if (is_non_ad_metal_name(name))
      return specific_atom_type::GenericMetal; //generic metal
    return specific_atom_type::NumTypes;
  } else {
    for (int i = 0; i < specific_atom_type::NumTypes; i++) {
      if (specific_atom_type::data[i].special_name == name) {
        return specific_atom_type::data[i].sm;
      }
    }
    return specific_atom_type::NumTypes;
  }
}

inline int num_atom_types() {
  return specific_atom_type::NumTypes;
}

inline bool is_hydrogen(smt t) {
  return t == specific_atom_type::Hydrogen || t == specific_atom_type::PolarHydrogen;
}

inline bool is_heteroatom(smt t) {
  assert(t < specific_atom_type::NumTypes);
  return specific_atom_type::data[t].ad_heteroatom;
}

inline double xs_radius(smt t) {
  assert(t < specific_atom_type::NumTypes);
  return specific_atom_type::data[t].xs_radius;
}

inline bool xs_is_hydrophobic(smt sm) {
  assert(sm < specific_atom_type::NumTypes);
  return specific_atom_type::data[sm].xs_hydrophobe;
}

inline bool xs_is_acceptor(smt sm) {
  assert(sm < specific_atom_type::NumTypes);
  return specific_atom_type::data[sm].xs_acceptor;
}

inline bool xs_is_donor(smt sm) {
  assert(sm < specific_atom_type::NumTypes);
  return specific_atom_type::data[sm].xs_donor;
}

inline bool xs_donor_acceptor(smt t1, smt t2) {
  return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(smt t1, smt t2) {
  return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}
