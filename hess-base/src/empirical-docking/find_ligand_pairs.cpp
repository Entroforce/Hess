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

#include "optimizable_molecule.h"
#include <queue>
using namespace std;

void calc_distances(vector<int> &d, hess::Molecule *lig, int start_id) {
  int n = lig->get_atoms_count();
  vector<int> v(n, 1);
  int min_id;
  int temp;
  int min;
  d[start_id] = 0;
  int iters = 0;
  do {
    iters++;
    min_id = 1e9;
    min = 1e9;
    for (int i = lig->vertexBegin(); i != lig->vertexEnd(); i = lig->vertexNext(i)) {
      if ((v[i] == 1) && (d[i] < min)) {
        min = d[i];
        min_id = i;
      }
    }
    if (min_id != 1e9) {
      ver atom_ver = lig->getVertex(min_id);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei_min = atom_ver.neiVertex(id);
        temp = min + 1;
        if (temp < d[nei_min]) {
          d[nei_min] = temp;
        }
      }
      v[min_id] = 0;
    }
    if (iters > 1e9)
      throw HessException("Error when going into the depth of a molecule when calculating distances to other atoms!");
  } while (min_id < 1e9);
}

void find_ligand_pairs(hess::Molecule *lig) {
  vector<int> processed;
  set<pair<int, int>> lock_pairs;
  vector<pair<int, int>> new_pairs;
  processed.clear();
  lock_pairs.clear();
  int n, cur = 0;
  n = lig->atoms.size();
  vector <int> was(n, -1);
  std::queue <int> q;
  for (int i = 0; i < n; i++) {
    if (was[i] != -1)
      continue;
    q.push(i);
    while (!q.empty()) {
      int v = q.front();
      q.pop();
      if (was[v] != -1)
        continue;
      was[v] = cur;
      ver atom_ver = lig->getVertex(v);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int nei = atom_ver.neiVertex(id);
        int bond_id = lig->findEdgeIndex(v, nei);
        if (!lig->get_bond(bond_id)->rotatable && was[nei] == -1)
          q.push(nei);
      }
    }
    cur++;
  }
  std::vector<int> neicomp;
  neicomp.resize(lig->vertexEnd());
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    if (lig->getAtomNumber(a_id) == 1)
      continue;
    vector<int> dists(n, 1e9);
    calc_distances(dists, lig, a_id);
    for (int b_id = lig->vertexNext(a_id); b_id != lig->vertexEnd(); b_id = lig->vertexNext(b_id)) {
      if (lig->getAtomNumber(b_id) == 1)
        continue;
      if (dists[b_id] >= 4 && was[a_id] != was[b_id]) {
        bool was_false_root = false;
        
        for (int anei = lig->getVertex(a_id).neiBegin(); anei != lig->getVertex(a_id).neiEnd(); anei = lig->getVertex(a_id).neiNext(anei)) {
          if (lig->get_bond(lig->getVertex(a_id).neiEdge(anei))->rotatable && was[lig->getVertex(a_id).neiVertex(anei)] == was[b_id]) {
            was_false_root = true;
            break;
          }
        }

        for (int bnei = lig->getVertex(b_id).neiBegin(); bnei != lig->getVertex(b_id).neiEnd(); bnei = lig->getVertex(b_id).neiNext(bnei)) {
          if (lig->get_bond(lig->getVertex(b_id).neiEdge(bnei))->rotatable && was[lig->getVertex(b_id).neiVertex(bnei)] == was[a_id]) {
            was_false_root = true;
            break;
          }
        }
      
        for (int i = 0; i < lig->vertexEnd(); i++)
          neicomp[i] = 0;

        for (int anei = lig->getVertex(a_id).neiBegin(); anei != lig->getVertex(a_id).neiEnd(); anei = lig->getVertex(a_id).neiNext(anei)) {
          if (lig->get_bond(lig->getVertex(a_id).neiEdge(anei))->rotatable)
            neicomp[was[lig->getVertex(a_id).neiVertex(anei)]] = 1;
        }

        for (int bnei = lig->getVertex(b_id).neiBegin(); bnei != lig->getVertex(b_id).neiEnd(); bnei = lig->getVertex(b_id).neiNext(bnei)) {
          if (lig->get_bond(lig->getVertex(b_id).neiEdge(bnei))->rotatable)
          {
            if (neicomp[was[lig->getVertex(b_id).neiVertex(bnei)]] == 1) {
              was_false_root = true;
              break;
            }
          }
        }
        if (was_false_root) {
          continue;
        }
        lock_pairs.insert(make_pair(a_id, b_id));
      }
    }
  }
      
  
  for (auto &p : lock_pairs) {
    lig->pairs.push_back(p);
  }
}
