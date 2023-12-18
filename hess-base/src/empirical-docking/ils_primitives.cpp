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
#include "transformation.h"
#include "LBFGS.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXf;
using namespace LBFGSpp;

bool check_conf(hess::Molecule *lig, simplified_tree& tr, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  transformation(lig, tr, encoding_inv, x);
  bool bull = !geometry_changed(lig);
  return bull;
}

bool check_exceeded_box_limits_start(hess::Molecule *lig, double size_x, double size_y, double size_z, simplified_tree& tr, int& ex_count, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  size_x /= 2;
  size_y /= 2;
  size_z /= 2;
  transformation(lig, tr, encoding_inv, x);
  ex_count = 0;
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    hess::Atom* a = lig->get_atom(a_id);
    double x = a->x;
    double y = a->y;
    double z = a->z;
    if (x > size_x || y > size_y || z > size_z) {
      ex_count += 1;
    } else if (x < -size_x || y < -size_y || z < -size_z) {
      ex_count += 1;
    }
  }
  if (ex_count >= 1)
    return true;
  return false;
}

bool check_exceeded_box_limits(hess::Molecule *lig, double size_x, double size_y, double size_z, simplified_tree& tr, int& ex_count, const vector<int>& encoding_inv, const Eigen::VectorXd& x) {
  size_x /= 2;
  size_y /= 2;
  size_z /= 2;
  transformation(lig, tr, encoding_inv, x);
  ex_count = 0;
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    hess::Atom* a = lig->get_atom(a_id);
    double x = a->x;
    double y = a->y;
    double z = a->z;
    if (x > size_x || y > size_y || z > size_z) {
      ex_count += 1;
    } else if (x < -size_x || y < -size_y || z < -size_z) {
      ex_count += 1;
    }
  }
  if (ex_count >= lig->get_atoms_count())
    //    if (ex_count >= 1)
    return true;
  return false;
}

double random_fl_ab(double min, double max) {
  double rand_num = (double) rand() / RAND_MAX * (max - min) + min;
  assert(rand_num >= min && min <= max);
  return rand_num;
}

int random_int_ab(int min, int max) {
  int rand_num = rand() % (max - min + 1) + min;
  assert(rand_num >= min && min <= max);
  return rand_num;
}

void random_in_sphere(double* norm_vec) {
  while (true) {
    double r1 = random_fl_ab(-1, 1);
    double r2 = random_fl_ab(-1, 1);
    double r3 = random_fl_ab(-1, 1);
    if (sqr(r1) + sqr(r2) + sqr(r3) < 1) {
      norm_vec[0] = r1;
      norm_vec[1] = r2;
      norm_vec[2] = r3;
      return;
    }
  }
}

void random_in_sphere(hess::Vec3d& vec) {
  while (true) {
    double r1 = random_fl_ab(-1, 1);
    double r2 = random_fl_ab(-1, 1);
    double r3 = random_fl_ab(-1, 1);
    if (sqr(r1) + sqr(r2) + sqr(r3) < 1) {
      vec[0] = r1;
      vec[1] = r2;
      vec[2] = r3;
      return;
    }
  }
}

void first_randomize_monte_carlo(Eigen::VectorXd& x, Optimizable_molecule& mol) {
  for (int i = 0; i < x.size() - 3; i++) {
    double ran_fl = random_fl_ab(-M_PI, M_PI);
    x[i] = ran_fl;
  }
  int corner_id = 0;
  for (int i = x.size() - 3; i < x.size(); i++) {
    double ran_fl = random_fl_ab(mol.corner1[corner_id], mol.corner2[corner_id]);
    x[i] = ran_fl;
    corner_id++;
  }
}

double gyration_radius(Eigen::VectorXd& x, Optimizable_molecule& mol) {
  hess::Molecule* lig = mol.ligand;
  transformation(lig, mol.tr, mol.encoding_inv, x);
  double acc = 0;
  unsigned counter = 0;
  int orign_id = mol.tr.rootid;
  hess::Atom* orign_atom = lig->get_atom(orign_id);
  const hess::Vec3d& orign_vec = orign_atom->get_vector();
  for (int a_id = lig->vertexBegin(); a_id != lig->vertexEnd(); a_id = lig->vertexNext(a_id)) {
    hess::Atom* a = lig->get_atom(a_id);
    if (lig->getAtomNumber(a_id) != 1) {
      const hess::Vec3d& a_vec = a->get_vector();
      acc += hess::Vec3d::dist_sqr(a_vec, orign_vec);
      counter++;
    }
  }
  return (counter > 0) ? std::sqrt(acc / counter) : 0;
}

void normalize_angle(double& x) { // subtract or add enough 2*M_PI's to make x be in [-M_PI, M_PI]
  if (x > 3 * M_PI) { // very large
    double n = (x - M_PI) / (2 * M_PI); // how many 2*M_PI's do you want to subtract?
    x -= 2 * M_PI * std::ceil(n); // ceil can be very slow, but this should not be called often
    normalize_angle(x);
  } else if (x < -3 * M_PI) { // very small
    double n = (-x - M_PI) / (2 * M_PI); // how many 2*M_PI's do you want to add?
    x += 2 * M_PI * std::ceil(n); // ceil can be very slow, but this should not be called often
    normalize_angle(x);
  } else if (x > M_PI) { // in (   M_PI, 3*M_PI]
    x -= 2 * M_PI;
  } else if (x < -M_PI) { // in [-3*M_PI,  -M_PI)
    x += 2 * M_PI;
  }
  assert(x >= -M_PI && x <= M_PI);
  // in [-M_PI, M_PI]
}

void angle_to_quaternion(hess::Vec3d& axis, double angle, double* qt) { // axis is assumed to be a unit vector
  normalize_angle(angle); // this is probably only necessary if angles can be very big
  double c = std::cos(angle / 2);
  double s = std::sin(angle / 2);
  qt[0] = c;
  qt[1] = s * axis[0];
  qt[2] = s * axis[1];
  qt[3] = s * axis[2];
}

void angle_to_quaternion(const hess::Vec3d& rotation, double* qt) {
  double angle = rotation.length();
  if (angle > epsilon_fl) {
    hess::Vec3d axis;
    axis.x = (1 / angle) * rotation.x;
    axis.y = (1 / angle) * rotation.y;
    axis.z = (1 / angle) * rotation.z;
    angle_to_quaternion(axis, angle, qt);
    return;
  }
  qt[0] = 1;
  qt[1] = 0;
  qt[2] = 0;
  qt[3] = 0;
}

void mutate_solve(Eigen::VectorXd& x, Optimizable_molecule& mol) { // 2A for position, similar amp for orientation, randomize torsion
  int rots_count = mol.rot_bonds_count;
  int solve_count = mol.solve_count;
  int entities = 2 + mol.rot_bonds_count;
  if (entities == 0) return;
  int which_int = random_int_ab(0, int(entities - 1));
  size_t which = size_t(which_int);
  if (which == 0) {
    double pos_change[3] = {0};
    random_in_sphere(pos_change);
    x[solve_count - 3] += pos_change[0] * mol.mut_amplitude;
    x[solve_count - 2] += pos_change[1] * mol.mut_amplitude;
    x[solve_count - 1] += pos_change[2] * mol.mut_amplitude;
    return;
  }
  --which;
  if (which == 0) {
    double g_rad = gyration_radius(x, mol);
    if (g_rad > epsilon_fl) {
      hess::Vec3d rotation;
      hess::Vec3d rotation_dis;
      random_in_sphere(rotation_dis);
      rotation[0] = mol.mut_amplitude / g_rad * rotation_dis[0];
      rotation[1] = mol.mut_amplitude / g_rad * rotation_dis[1];
      rotation[2] = mol.mut_amplitude / g_rad * rotation_dis[2];
      double qt[4] = {0};
      angle_to_quaternion(rotation, qt);

      double q0 = qt[0];
      double q1 = qt[1];
      double q2 = qt[2];
      double q3 = qt[3];
      double rotx_dur = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
      double roty_dur = asin(2 * (q0 * q2 - q3 * q1));
      double rotz_dur = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));

      double curr_rotx = x[solve_count - 6];
      double curr_roty = x[solve_count - 5];
      double curr_rotz = x[solve_count - 4];
      Matrix4 M1 = global_rotate_by_x(curr_rotx);
      Matrix4 M2 = global_rotate_by_y(curr_roty);
      Matrix4 M3 = global_rotate_by_z(curr_rotz);
      Matrix4 M_curr = mlt(M3, M2);
      M_curr = mlt(M_curr, M1);
      M1 = global_rotate_by_x(rotx_dur);
      M2 = global_rotate_by_y(roty_dur);
      M3 = global_rotate_by_z(rotz_dur);
      Matrix4 M_dur = mlt(M3, M2);
      M_dur = mlt(M_dur, M1);
      Matrix4 M_new = mlt(M_curr, M_dur);
      double angles[3] = {0};
      rotate_matrix_to_euler_ZYX(M_new, angles);
      double x_angle = angles[0];
      double y_angle = angles[1];
      double z_angle = angles[2];
      double new_rotx = x_angle;
      double new_roty = y_angle;
      double new_rotz = z_angle;
      normalize_angle(new_rotx);
      normalize_angle(new_roty);
      normalize_angle(new_rotz);
      x[solve_count - 6] = new_rotx;
      x[solve_count - 5] = new_roty;
      x[solve_count - 4] = new_rotz;
    }
    return;
  }
  --which;
  if (which < rots_count) {
    x[which] = random_fl_ab(-M_PI, M_PI);
    return;
  }
  assert(false);
}

bool metropolis_accept(double old_f, double new_f, double temperature) {
  if (new_f < old_f) return true;
  const double acceptance_probability = std::exp((old_f - new_f) / temperature);
  return random_fl_ab(0, 1) < acceptance_probability;
}

double rmsd_upper_bound(const vector<hess::Vec3d>& a, const vector<hess::Vec3d>& b) {
  double acc = 0;
  for (int i = 0; i < a.size(); i++) {
    const hess::Vec3d& a_vec = a[i];
    const hess::Vec3d& b_vec = b[i];
    acc += hess::Vec3d::dist_sqr(a[i], b[i]);
  }
  return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}

std::pair<size_t, double> find_closest(const vector<hess::Vec3d>& a, const vector<mc_out>& b) {
  std::pair<size_t, double> tmp(b.size(), max_fl);
  for (int i = 0; i < b.size(); i++) {
    const vector<hess::Vec3d>& b_coords = b[i].get_coords();
    double res = rmsd_upper_bound(a, b_coords);
    if (i == 0 || res < tmp.second)
      tmp = std::pair<size_t, double>(i, res);
  }
  return tmp;
}

void add_to_output_container(vector<mc_out>& outs, VectorXd& tmp, Optimizable_molecule& mol, double e) {
  transformation(mol.ligand, mol.tr, mol.encoding_inv, tmp);
  mol.ligand->calc_coords();
  const vector<hess::Vec3d>& coords = mol.ligand->coords;
  std::pair<size_t, double> closest_rmsd = find_closest(coords, outs);
  if (closest_rmsd.first < outs.size() && closest_rmsd.second < mol.min_rmsd) { // have a very similar one
    if (e < outs[closest_rmsd.first].get_energy()) { // the new one is better, apparently
      outs[closest_rmsd.first] = mc_out(coords, tmp, e);
    }
  } else { // nothing similar
    if (outs.size() < mol.num_saved_mins)
      outs.push_back(mc_out(coords, tmp, e)); // the last one had the worst energy - replacing 
    else
      if (!outs.empty() && e < outs.back().get_energy())
      outs.back() = mc_out(coords, tmp, e);
  }
  sort(outs.begin(), outs.end(), [ ](const auto& l, const auto& r) {
    return l.e < r.e; });
}

void add_to_output_container_for_res(vector<mc_out>& outs, const VectorXd& tmp, Optimizable_molecule& mol, double e, const vector<hess::Vec3d>& coords) {
  std::pair<size_t, double> closest_rmsd = find_closest(coords, outs);
  if (closest_rmsd.first < outs.size() && closest_rmsd.second < mol.min_rmsd) { // have a very similar one
    if (e < outs[closest_rmsd.first].get_energy()) { // the new one is better, apparently
      outs[closest_rmsd.first] = mc_out(coords, tmp, e);
    }
  } else { // nothing similar
    if (outs.size() < mol.num_saved_mins)
      outs.push_back(mc_out(coords, tmp, e)); // the last one had the worst energy - replacing 
    else
      if (!outs.empty() && e < outs.back().get_energy())
      outs.back() = mc_out(coords, tmp, e);
  }
  sort(outs.begin(), outs.end(), [ ](const auto& l, const auto& r) {
    return l.e < r.e; });
}

void merge_output_containers(vector<mc_out>& all_outs, vector<mc_out>& outs, Optimizable_molecule& mol) {
  for (int i = 0; i < outs.size(); i++) {
    const VectorXd& tmp = outs[i].get_solve();
    const vector<hess::Vec3d>& coords = outs[i].get_coords();
    double tmp_e = outs[i].get_energy();
    add_to_output_container_for_res(all_outs, tmp, mol, tmp_e, coords); // 20 - max size
  }
}