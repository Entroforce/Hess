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

void move_box(vector<double>& box, double x_diff, double y_diff, double z_diff) {
  box[0] = x_diff;
  box[1] = y_diff;
  box[2] = z_diff;
}

void calc_global_shifts(vector<vector<double>>&global_shifts, Optimizable_molecule& mol, int niter) {
  double size_x = mol.size_x;
  double size_y = mol.size_y;
  double size_z = mol.size_z;
  double lig_box_size_x = mol.ligand_box[0] + mol.ligand_box[1];
  double lig_box_size_y = mol.ligand_box[2] + mol.ligand_box[3];
  double lig_box_size_z = mol.ligand_box[4] + mol.ligand_box[5];
  double lx = mol.ligand_box[0];
  double rx = mol.ligand_box[1];
  double ly = mol.ligand_box[2];
  double ry = mol.ligand_box[3];
  double lz = mol.ligand_box[4];
  double rz = mol.ligand_box[5];
  if (lig_box_size_x >= size_x) {
    if (lx >= size_x / 2) {
      lx = size_x / 4;
    }
    if (rx >= size_x / 2) {
      rx = size_x / 4;
    }
  }
  if (lig_box_size_y >= size_y) {
    if (ly >= size_y / 2) {
      ly = size_y / 4;
    }
    if (ry >= size_y / 2) {
      ry = size_y / 4;
    }
  }
  if (lig_box_size_z >= size_z) {
    if (lz >= size_z / 2) {
      lz = size_z / 4;
    }
    if (rz >= size_z / 2) {
      rz = size_z / 4;
    }
  }
  lig_box_size_x = lx + rx;
  lig_box_size_y = ly + ry;
  lig_box_size_z = lz + rz;
  int STEPS_COUNT = niter;
  double axis_steps = std::pow(STEPS_COUNT, 1.0 / 3.0);
  vector<double> lig_center(3);
  lig_center[0] = mol.ligand_center[0];
  lig_center[1] = mol.ligand_center[1];
  lig_center[2] = mol.ligand_center[2];
  double xd = (size_x - lig_box_size_x) / axis_steps;
  double yd = (size_y - lig_box_size_y) / axis_steps;
  double zd = (size_z - lig_box_size_z) / axis_steps;
  double xd_i = -size_x / 2 + lx;
  double yd_i = -size_y / 2 + ly;
  double zd_i = -size_z / 2 + lz;
  int z_c = 0, y_c = 0, x_c = 0;
  while (xd_i < size_x / 2 - rx) {
    yd_i = -size_y / 2 + ly;
    while (yd_i < size_y / 2 - ry) {
      zd_i = -size_z / 2 + lz;
      while (zd_i < size_z / 2 - rz) {
        move_box(lig_center, xd_i, yd_i, zd_i);
        global_shifts.push_back(lig_center);
        zd_i += (zd);
        z_c += 1;
      }
      yd_i += (yd);
      y_c += 1;
    }
    xd_i += (xd);
    x_c += 1;
  }
}

void moveLigandToBoxCenter(double xc, double yc, double zc, void* bestmol_v) {
  hess::Molecule* bestmol = (hess::Molecule*)bestmol_v;
  for (int a_id = bestmol->vertexBegin(); a_id != bestmol->vertexEnd(); a_id = bestmol->vertexNext(a_id)) {
    hess::Atom* a = bestmol->get_atom(a_id);
    a->x += xc;
    a->y += yc;
    a->z += zc;
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
