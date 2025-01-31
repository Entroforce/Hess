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

#include "vector3.h"
#include <ctime>

using namespace indigo;

bool equal(const hess::Vec3d& a, const hess::Vec3d& b) {
  return ((a.x == b.x) &&
          (a.y == b.y) &&
          (a.z == b.z));
}

bool is_zero(const hess::Vec3d& v) {
  return ((v.x == 0.0) &&
          (v.y == 0.0) &&
          (v.z == 0.0));
}

void random_vector(hess::Vec3d& a) {
  double l;
  srand(static_cast<unsigned int> (time(nullptr)));
  double r1 = double(rand()) / RAND_MAX;
  double r2 = double(rand()) / RAND_MAX;
  double r3 = double(rand()) / RAND_MAX;
  do {
    a.set(r1 - 0.5, r2 - 0.5, r3 - 0.5);
    l = a.lengthSqr();
  } while ((l > 1.0) || (l < 1e-4));
  a.normalize();
}
