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

#include <sys/stat.h>
#include "bond/bond_order_determination.h"
#include "parser/parser.h"
#include "scoring_function/features.h"
#include "scoring_function/scoring.h"

#include "array3d.h"
#include "fast_gradient.h"
#include "ils_primitives.h"
#include "global_optimizer.h"
#include "hess.h"
#include "vector3.h"
#include "logger.h"
#include "base_cpp/output.h"
#include "molecule/molfile_saver.h"

using namespace std;
using namespace hess;

thread_local static std::string _error_message;

struct HessObject {

  HessObject(Molecule* _mol = nullptr, Optimizable_molecule* _opt_mol = nullptr, Parser* _parser = nullptr) : mol(_mol), opt_mol(_opt_mol), parser(_parser) {
  }

  ~HessObject() {
    if (mol != nullptr) {
      delete mol;
    } else if (opt_mol != nullptr) {
      delete opt_mol;
    } else if (parser != nullptr) {
      delete parser;
    }
  }

  Molecule* mol;
  Optimizable_molecule* opt_mol;
  Parser* parser;
};

Molecule* to_molecule(void* object) {
  HessObject* ho = (HessObject*) object;
  if (ho->mol == nullptr) {
    throw new HessException("The object is not a molecule!\n");
  }
  return ho->mol;
}

Optimizable_molecule* to_optmolecule(void* object) {
  HessObject* ho = (HessObject*) object;
  if (ho->opt_mol == nullptr) {
    throw new HessException("The object is not a optimizable molecule.\n");
  }
  return ho->opt_mol;
}

Parser* to_parser(void* object) {
  HessObject* ho = (HessObject*) object;
  if (ho->parser == nullptr) {
    throw new HessException("The object is not a parser.\n");
  }
  return ho->parser;
}

void hessSetLogFile(const char * path) {
  Logger::get_instance().set_path(path);
}

FILE* hessGetStream() {
  return Logger::get_instance().get_stream();
}

void hessSetStream(FILE* stream) {
  Logger::get_instance().set_stream(stream);
}

void hessInit() {
  indigoAllocSessionId();
}

void *hessCreateParser() {
  return new HessObject(nullptr, nullptr, new hess::Parser());
}

void* hessLoadMolecule(void *parser, const char *file_name) {
  try {
    return new HessObject(to_parser(parser)->parse_molecule(file_name), nullptr, nullptr);
  } catch (HessException &e) {
    _error_message = e.what();
    return NULL;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return NULL;
  } catch (exception &e) {
    _error_message = e.what();
    return NULL;
  }
}

int hessSaveMol(void *mol, const char *path) {
  try {
    indigo::FileOutput output(path);
    indigo::MolfileSaver saver(output);
    saver.saveBaseMolecule(*to_molecule(mol));
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

int hessDeleteFormalCharges(void *mol) {
  hess::Molecule *molecule = to_molecule(mol);
  try {
    for (int a_id = molecule->vertexBegin(); a_id != molecule->vertexEnd(); a_id = molecule->vertexNext(a_id)) {
      molecule->setAtomCharge(a_id, 0);
    }
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }

}

int hessSaveSdf(void* opt_molecule, const char *path) {
  try {
    Optimizable_molecule *opt = to_optmolecule(opt_molecule);
    vector<hess::Molecule*> res_mols = opt->result_mols;
    indigo::FileOutput output(path);
    indigo::MolfileSaver saver(output);
    saver.mode = indigo::MolfileSaver::MODE_3000;
    for (int i = 0; i < res_mols.size(); i++) {
      hess::Molecule* mol = res_mols[i];
      saver.saveBaseMolecule(*mol);
      output.printfCR("$$$$");
    }
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

void hessSetPh(void *mol, double pH) {
  hess::Molecule *molecule = to_molecule(mol);
  molecule->pH = pH;
}

int hessDeleteHydrogens(void *mol) {
  hess::Molecule *molecule = to_molecule(mol);
  try {

    vector<int> to_remove;

    for (int i = molecule->vertexBegin(); i != molecule->vertexEnd(); i = molecule->vertexNext(i))
      if (molecule->getAtomNumber(i) == 1)
        to_remove.push_back(i);

    for (int i = 0; i < to_remove.size(); i++)
      molecule->remove_atom(to_remove[i]);
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

int hessAssignBonds(void* old_mol) {
  hess::Molecule *mol = to_molecule(old_mol);
  try {
    determine_connectivity(mol);
    delete_bonds_exceeding_valence(mol);
    hydridize_atoms(mol);
    vector<vector<int>> rings_bonds;
    vector<vector<int>> rings_atoms;
    extract_rings(rings_bonds, rings_atoms, mol);
    mol->set_rings(rings_atoms);
    check_torsion_angles(mol, rings_atoms);
    antialiasing(mol);
    functional_groups(mol);
    for (int i = 0; i < rings_bonds.size(); i++)
      for (int j = 0; j < rings_atoms[i].size(); j++)
        mol->get_atom(rings_atoms[i][j])->in_ring = true;

    assign_bond_orders(mol, false);
    determine_bond_orders_in_aromatic_rings(mol, rings_bonds, rings_atoms);
    mol->set_aromatic_rings(move(rings_atoms));
    assign_bond_orders(mol);
    hybridize(mol);
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

int hessProtonate(void* mol) {
  try {
    hess::Molecule *molecule = to_molecule(mol);
    if (molecule->pH >= 0) {
      transform_Ph(molecule);
    }
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

int hessAddHydrogens(void *mol) {
  try {
    hess::Molecule *atoms = to_molecule(mol);
    // add new implicit hydrogens
    unfold_hydrogens(atoms);
    atoms->charge_calculate();
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

void hessCalcAutobox(void* ligand_v, double* box) {
  hess::Molecule *ligand = to_molecule(ligand_v);
  double size_x, size_y, size_z, xc, yc, zc;
  double min_x = 1e9, min_y = 1e9, min_z = 1e9;
  double max_x = -1e9, max_y = -1e9, max_z = -1e9;
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    hess::Atom *a = ligand->get_atom(a_id);
    double x = a->x;
    double y = a->y;
    double z = a->z;
    if (x < min_x)
      min_x = x;
    if (x > max_x)
      max_x = x;
    if (y < min_y)
      min_y = y;
    if (y > max_y)
      max_y = y;
    if (z < min_z)
      min_z = z;
    if (z > max_z)
      max_z = z;
  }
  xc = (min_x + max_x) / 2;
  yc = (min_y + max_y) / 2;
  zc = (min_z + max_z) / 2;
  size_x = fabs(max_x - min_x) + 8;
  size_y = fabs(max_y - min_y) + 8;
  size_z = fabs(max_z - min_z) + 8;
  box[0] = xc;
  box[1] = yc;
  box[2] = zc;
  box[3] = size_x;
  box[4] = size_y;
  box[5] = size_z;
}

void hessWriteScoreOnly(void* opt_mol) {
  Optimizable_molecule* opt = to_optmolecule(opt_mol);
  hess::Molecule* lig = opt->ligand;
  hess::Molecule* rec = opt->receptor;
  double e = calc_affinity(lig, rec, opt->scoring);
  double e_intra = calc_intramolecular_energy(lig, rec, opt->scoring);
  fprintf(hessGetStream(), "Affinity: %7.3f (kcal/mol)\n", e);
  fprintf(hessGetStream(), "Intramolecular energy: %7.3f (kcal/mol)\n", e_intra);

}

void* hessMakeOptimizableMolecule(void* ligand_v, void* rec_v, double* box, const char* optimize, double granularity, unsigned seed) {
  try {
    return new HessObject(nullptr, new Optimizable_molecule(to_molecule(ligand_v), to_molecule(rec_v), box, optimize, granularity, seed), nullptr);
  } catch (HessException &e) {
    _error_message = e.what();
    return NULL;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return NULL;
  } catch (exception &e) {
    _error_message = e.what();
    return NULL;
  }
}

void hessFillGrid(void* opt_v, int grid_deriv_flag) {
  Optimizable_molecule *opt = to_optmolecule(opt_v);
  if (!grid_deriv_flag) {
    opt->fill_grid();
  } else {
    opt->fill_grid_with_derivs();
    opt->grid_with_derivs = true;
  }
}

void hessFillGridWithDerivs(void* opt_v) {
  Optimizable_molecule* opt = (Optimizable_molecule*) opt_v;
  opt->fill_grid_with_derivs();
  opt->grid_with_derivs = true;
}

int hessRunOptimize(int number_of_iterations, int depth, void* opt_v, double* result_array, int tops_count) {
  Optimizable_molecule *opt = to_optmolecule(opt_v);
  int error_flag = 0;
  if (opt->optimize == "mc") {
    error_flag = hessRunRandomIls(number_of_iterations, depth, opt_v, result_array, tops_count);
  } else if (opt->optimize == "mc_metropolis") {
    error_flag = hessRunMonteCarlo(number_of_iterations, opt_v, result_array, tops_count);
  }
  return error_flag;
}

int hessRunRandomIls(int number_of_iterations, int depth, void* opt_v, double* result_array, int tops_count) {
  try {
    Optimizable_molecule *opt = to_optmolecule(opt_v);
    hess::Molecule* ord_lig = opt->ligand;
    hess::Molecule* ord_rec = opt->receptor;
    const char* scoring_type = opt->scoring;
    double dif = 0.1;
    unsigned max_bfgs_iterations = unsigned((25.0 + ord_lig->get_atoms_count()) / 3.0);
    opt->set_max_bfgs_iterations(max_bfgs_iterations);
    vector< pair< Eigen::VectorXd, pair<double, double> > >result_pairs;
    srand(opt->seed);
    fprintf(hessGetStream(), "Seed: %d\n", opt->seed);
    fprintf(hessGetStream(), "Starting Monte Carlo. Number of iterations: %d. Depth: %d\n", number_of_iterations, depth);
    for (int i = 0; i < number_of_iterations; i++) {
      Eigen::VectorXd ils_result;

      ils_result = ils_random(*opt, depth, dif);

      if (!check_conf(ord_lig, opt->tr, opt->encoding_inv, ils_result)) {
        fprintf(stderr, "Form changed.\n");
      }
      pair<double, double> res = calc_energy_for_result(ord_lig, ord_rec, opt->tr, opt->in, opt->encoding_inv, ils_result, scoring_type);
      fprintf(hessGetStream(), "Iteration %3i Inter energy: %7.3f Intra energy: %7.3f Sum: %7.3f\n",
              i + 1, res.first, res.second, res.first + res.second);
      result_pairs.push_back({ils_result,
        {res.first, res.second}});
    }
    sort_configurations(result_pairs);
    opt->tops_count = min((unsigned) result_pairs.size(), (unsigned) tops_count);
    form_ils_results(result_pairs, opt);
    result_array[0] = result_pairs[0].second.first;
    result_array[1] = result_pairs[0].second.second;
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

int hessRunMonteCarlo(int number_of_iterations, void* opt_v, double* result_array, int tops_count) {
  try {
    Optimizable_molecule *opt = to_optmolecule(opt_v);
    const char* scoring_type = opt->scoring;
    hess::Molecule* ord_lig = opt->ligand;
    hess::Molecule* ord_rec = opt->receptor;
    unsigned max_bfgs_iterations = unsigned((25.0 + ord_lig->get_atoms_count()) / 3.0);
    opt->set_max_bfgs_iterations(max_bfgs_iterations);
    vector< pair< Eigen::VectorXd, pair<double, double> > >result_pairs;
    vector<mc_out> mc_all_results;
    srand(opt->seed);
    fprintf(hessGetStream(), "Seed: %d\n", opt->seed);
    fprintf(hessGetStream(), "Starting Monte Carlo with Metropolis. Number of steps: %d. Number of iterations: %d\n", opt->num_steps, number_of_iterations);
    for (int i = 0; i < number_of_iterations; i++) {
      vector<mc_out> mc_results;
      monte_carlo(*opt, mc_results);
      fprintf(hessGetStream(), "Iteration %3i Current top energy: %7.3g\n", i + 1, mc_results[0].get_energy());
      merge_output_containers(mc_all_results, mc_results, *opt);
    }
    for (int i = 0; i < opt->num_saved_mins; i++) {
      Eigen::VectorXd mc_solve = VectorXd::Zero(opt->encoding.size());
      mc_solve = mc_all_results[i].solve;
      pair<double, double> res = calc_energy_for_result(ord_lig, ord_rec, opt->tr, opt->in, opt->encoding_inv, mc_solve, scoring_type);
      result_pairs.push_back({mc_solve,
        {res.first, res.second}});
    }
    sort_configurations(result_pairs);
    opt->tops_count = min((unsigned) result_pairs.size(), (unsigned) tops_count);
    form_ils_results(result_pairs, opt);
    result_array[0] = result_pairs[0].second.first;
    result_array[1] = result_pairs[0].second.second;
    return 0;
  } catch (HessException &e) {
    _error_message = e.what();
    return -1;
  } catch (indigo::Exception &e) {
    _error_message = e.what();
    return -1;
  } catch (exception &e) {
    _error_message = e.what();
    return -1;
  }
}

void hessDestroy(void *object) {
  delete (HessObject*) object;
}

const char *hessGetLastError() {
  return _error_message.c_str();
}
