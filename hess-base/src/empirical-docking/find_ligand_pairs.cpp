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
using namespace std;

void calc_distances(vector<int> &d, hess::Molecule *lig, int start_id, vector <int>& rot_path) {
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
  //find path for rotatable flag
  vector<int> processed;
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    if (a_id != start_id)
      processed.push_back(a_id);
  }
  for (int i : processed) {
    int end = i;
    int w = d[end];
    bool was_rot = false;
    iters = 0;
    while (end != start_id) {
      iters++;
      for (int id = lig->getVertex(end).neiBegin(); id != lig->getVertex(end).neiEnd();
               id = lig->getVertex(end).neiNext(id)) {
        int nei_id = lig->getVertex(end).neiVertex(id);
        int temp = w - 1;
        if (temp == d[nei_id]) {
          if (lig->get_bond(lig->findEdgeIndex(end, nei_id))->rotatable) {
            was_rot = true;
            break;
          }
          w = temp;
          end = nei_id;
          break;
        }
      }
      if (was_rot)
        break;
      if (iters > 1e9)
        throw HessException("Error when going into the depth of a molecule when calculating distances to other atoms!");
    }
    if (was_rot) {
      rot_path[i] = 1;
    }
  }
}


// if a_id and b_id relate to same ring (m.b. b_id connect to ring)

bool in_same_rings(hess::Molecule *lig, int a_id, int b_id) {
  int a_ring_id = -1;
  vector <vector<int>> rings = lig->get_rings();
  hess::Atom* a = lig->get_atom(a_id);
  if (lig->getVertex(a_id).degree() == 1 && !a->in_ring) {
    ver atom_ver = lig->getVertex(a_id);
    for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
      int nei_id = atom_ver.neiVertex(id);
      hess::Atom* nei = lig->get_atom(nei_id);
      if (nei->in_ring) {
        int ring_id = 0;
        for (vector<int> &ring : rings) {
          for (int ring_atom_id : ring) {
            if (ring_atom_id == nei_id) {
              a_ring_id = ring_id;
              break;
            }
          }
          if (a_ring_id != -1)
            break;
          ring_id++;
        }
      }
    }
  }
  if (!a->in_ring)
    return false;
  if (a_ring_id == -1) {
    int ring_id = 0;
    for (vector<int> &ring : rings) {
      for (int ring_atom_id : ring) {
        if (ring_atom_id == a_id) {
          a_ring_id = ring_id;
          break;
        }
      }
      if (a_ring_id != -1)
        break;
      ring_id++;
    }
  }
  assert(a_ring_id != -1);
  vector<int> a_ring = rings[a_ring_id];
  hess::Atom* b = lig->get_atom(b_id);
  if (b->in_ring) {
    for (int ring_id : a_ring) {
      if (ring_id == b_id) {
        return true;
      }
    }
  }
  ver atom_ver = lig->getVertex(b_id);
  for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
    int nei_id = atom_ver.neiVertex(id);
    for (int ring_id : a_ring) {
      if (ring_id == nei_id) {
        return true;
      }
    }
  }
  return false;
}

void find_ligand_pairs(hess::Molecule *lig) {
  int count = 0;
  vector<int> processed;
  set<pair<int, int>> lock_pairs;
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    hess::Atom* a = lig->get_atom(a_id);
    if (lig->getAtomNumber(a_id) == 1)
      continue;
    processed.push_back(a_id);
  }
  int n = lig->get_atoms_count();
  for (int p : processed) {
    vector<int> visited(n, 0);
    vector<int> dists(n, 1e9);
    vector <int> rot_path(n, 0);
    std::stack<int> s;
    s.push(p);
    calc_distances(dists, lig, p, rot_path);
    int iters = 0;
    while (!s.empty()) {
      if (iters > 1e9)
        throw HessException("Error when going into the depth of a molecule when searching for pairs of atoms in a ligand!");
      int u = s.top();
      visited[u] = 1;
      s.pop();
      ver atom_ver = lig->getVertex(u);
      for (int id = atom_ver.neiBegin(); id != atom_ver.neiEnd(); id = atom_ver.neiNext(id)) {
        int v = atom_ver.neiVertex(id);
        if (lig->getAtomNumber(v) == 1)
          continue;
        if (visited[v]) {
          continue;
        }
        if (dists[v] >= 4 && !in_same_rings(lig, p, v) && !in_same_rings(lig, v, p) && lock_pairs.find(make_pair(v, p)) == lock_pairs.end()
                && lock_pairs.find(make_pair(p, v)) == lock_pairs.end() && rot_path[v]) {
          count++;
          lock_pairs.insert(make_pair(p, v));
        }
        s.push(v);
      }
      iters++;
    }
  }
  for (auto &p : lock_pairs)
    lig->pairs.push_back(p);
}
