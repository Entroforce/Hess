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

#include "global_optimizer.h"
#include "logger.h"
#include "transformation.h"

using namespace std;

void sort_configurations(vector<pair<Eigen::VectorXd, pair<double, double>>>& result_pairs) {
  sort(result_pairs.begin(), result_pairs.end(), [ ](const auto& l, const auto& r) {
    return l.second.first + l.second.second < r.second.first + r.second.second; });
}

vector<double> calc_autobox(hess::Molecule *ligand) {
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
  printf("xc = %lf, yc = %lf, zc = %lf, sizex = %lf, sizey = %lf, sizez = %lf\n", xc, yc, zc, size_x, size_y, size_z);
  return {xc, yc, zc, size_x, size_y, size_z};
}

void calc_center_mass_frag(hess::Molecule *pd, simplified_tree &tr, int v, int& count, hess::Vec3d& average) {
  average.x += pd->atoms[v].x;
  average.y += pd->atoms[v].y;
  average.z += pd->atoms[v].z;
  count += 1;
  for (auto child : tr.children[v]) {
    if (!tr.bond_to_parent[child]) {
      calc_center_mass_frag(pd, tr, child, count, average);
    }
  }
}

void move_ligand_to_center(double* ligand_center, hess::Molecule *ligand) {
  for (int a_id = ligand->vertexBegin(); a_id != ligand->vertexEnd(); a_id = ligand->vertexNext(a_id)) {
    hess::Atom* a = ligand->get_atom(a_id);
    a->x -= ligand_center[0];
    a->y -= ligand_center[1];
    a->z -= ligand_center[2];
    a->x_orign -= ligand_center[0];
    a->y_orign -= ligand_center[1];
    a->z_orign -= ligand_center[2];
  }
}

void moveLigandToBoxCenter(double xc, double yc, double zc, void* bestmol_v) {
  hess::Molecule* bestmol = (hess::Molecule*)bestmol_v;
  for (int a_id = bestmol->vertexBegin(); a_id != bestmol->vertexEnd(); a_id = bestmol->vertexNext(a_id)) {
    hess::Atom* a = bestmol->get_atom(a_id);
    a->x += xc;
    a->y += yc;
    a->z += zc;
    bestmol->setAtomXyz(a_id, a->x, a->y, a->z);
  }
}

void form_ils_results(vector<pair<Eigen::VectorXd, pair<double, double>>>& result_pairs_sort, Optimizable_molecule* opt) {
  if (!result_pairs_sort.size()) {
    throw HessException("All optimal conformations of the ligand changed their shape. Try more iterations");
  }
  double top_intra = result_pairs_sort[0].second.second;
  fprintf(hess::Logger::get_instance().get_stream(), "\nResults:\n");
  for (int i = 0; i < opt->tops_count; i++) {
    fprintf(hess::Logger::get_instance().get_stream(), "Inter energy: %7.3f Intra energy: %7.3f Sum: %7.3f Sum-top intra: %7.3f\n",
            result_pairs_sort[i].second.first, result_pairs_sort[i].second.second,
            result_pairs_sort[i].second.first + result_pairs_sort[i].second.second,
            result_pairs_sort[i].second.first + result_pairs_sort[i].second.second - top_intra);
    hess::Molecule* opt_mol = new hess::Molecule(*((hess::Molecule*)opt->ligand));
    transformation(opt_mol, opt->tr, opt->encoding_inv, result_pairs_sort[i].first);
    moveLigandToBoxCenter(opt->xc, opt->yc, opt->zc, opt_mol);
    opt->result_mols.push_back(opt_mol);
  }
  fprintf(hess::Logger::get_instance().get_stream(), "\n");
}

void dfs_insubtree(int base, vector<vector<int>>&insubt, hess::Molecule *pd, simplified_tree &tr, int v, int tf) {
  int tf1 = tf;
  if (v == base)
    tf1 = 1;
  insubt[base][v] = tf1;
  for (auto child : tr.children[v])
    dfs_insubtree(base, insubt, pd, tr, child, tf1);
}

void set_insubtree(hess::Molecule* lig) {
  int root_id = lig->vertexBegin();
  simplified_tree tr(lig->get_atoms_count());
  tree_init(tr, lig);
  vector<vector<int>> insubt(tr.sz, vector<int>(tr.sz));
  for (int i = 0; i < tr.sz; i++)
    dfs_insubtree(i, insubt, lig, tr, tr.rootid, 0);
  lig->set_insubt(move(insubt));
}

void set_parents_to_tree(hess::Molecule* lig) {
  simplified_tree tr(lig->get_atoms_count());
  tree_init(tr, lig);
  vector<int> rotsmap;
  for (int i = 0; i < tr.bond_to_parent.size(); i++) {
    if (tr.bond_to_parent[i]) {
      rotsmap.push_back(i);
    }
  }
  vector<int> tree_parents(rotsmap.size());
  for (int i = 0; i < rotsmap.size(); i++) {
    int pos = rotsmap[i];
    int father = -1;
    for (int i1 = 0; i1 < tr.sz; i1++) {
      for (int j : tr.children[i1]) {
        if (j == pos) {
          father = i1;
        }
      }
    }
    tree_parents[i] = father;
  }
  lig->set_tree_parents(move(tree_parents));
}

pair<double, double> calc_energy_for_result(hess::Molecule* lig, hess::Molecule* prot, simplified_tree &tr, ConfIndependentInputs& in, const vector<int>& encoding_inv, const Eigen::VectorXd& x, const char* scoring) {
  transformation(lig, tr, encoding_inv, x);
  double inter_e = calc_affinity_with_confind(lig, prot, in, scoring);
  double intra_e = calc_intramolecular_energy(lig, prot, scoring);
  return {inter_e, intra_e};
}

double calc_energy(hess::Molecule *lig, hess::Molecule *prod, ConfIndependentInputs& in, const char* scoring) {
  return calc_affinity_with_confind(lig, prod, in, scoring);
}
