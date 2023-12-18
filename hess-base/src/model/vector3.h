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

#include <math.h>
#include <random>
#include "math/algebra.h"
#include "assert.h"

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180.0/M_PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI/180.0)
#endif

inline double sqr(double x) {
  return x * x;
}

namespace hess {

  struct Vec3d {

    Vec3d() : x(0), y(0), z(0) {
    }

    Vec3d(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {
    }

    Vec3d(double arr[3]) : x(arr[0]), y(arr[1]), z(arr[2]) {
    }

    Vec3d(const Vec3d& v) : x(v.x), y(v.y), z(v.z) {
    }

    Vec3d(Vec3d&& v) : x(v.x), y(v.y), z(v.z) {
    }

    Vec3d& operator=(const Vec3d& other) {
      x = other.x;
      y = other.y;
      z = other.z;
      return *this;
    }

    Vec3d& operator-=(const Vec3d& other) {
      x -= other.x;
      y -= other.y;
      z -= other.z;
      return *this;
    }

    double& operator[](int i) {
      assert(i < 3);
      if (i == 0)
        return this->x;
      if (i == 1)
        return this->y;
      return this->z;
    }

    Vec3d& operator=(Vec3d&& other) {
      x = other.x;
      y = other.y;
      z = other.z;
      return *this;
    }

    double x, y, z;

    static double dot(const Vec3d& a, const Vec3d& b) {
      return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    static double dist_sqr(const Vec3d& a, const Vec3d& b) {
      return sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z);
    }

    void elementwise_product(const Vec3d& a, const Vec3d& b) {
      x = a.x * b.x;
      y = a.y * b.y;
      z = a.z * b.z;

    }

    bool normalize() {
      double l = lengthSqr();

      if (l < indigo::EPSILON * indigo::EPSILON)
        return false;

      l = (double) sqrt(l);

      x /= l;
      y /= l;
      z /= l;

      return true;
    }

    void set(double xx, double yy, double zz) {
      x = xx;
      y = yy;
      z = zz;
    }

    void copy(const Vec3d& a) {
      x = a.x;
      y = a.y;
      z = a.z;
    }

    void zero() {
      x = 0;
      y = 0;
      z = 0;
    }

    void clear() {
      zero();
    }

    void negate() {
      x = -x;
      y = -y;
      z = -z;
    }

    void negation(const Vec3d& v) {
      x = -v.x;
      y = -v.y;
      z = -v.z;
    }

    void add(const Vec3d& v) {
      x += v.x;
      y += v.y;
      z += v.z;
    }

    void sum(const Vec3d& a, const Vec3d& b) {
      x = a.x + b.x;
      y = a.y + b.y;
      z = a.z + b.z;
    }

    void sub(const Vec3d& v) {
      x -= v.x;
      y -= v.y;
      z -= v.z;
    }

    void diff(const Vec3d& a, const Vec3d& b) {
      x = a.x - b.x;
      y = a.y - b.y;
      z = a.z - b.z;
    }

    void min(const Vec3d& a) {
      x = std::min(x, a.x);
      y = std::min(y, a.y);
      z = std::min(z, a.z);
    }

    void max(const Vec3d& a) {
      x = std::max(x, a.x);
      y = std::max(y, a.y);
      z = std::max(z, a.z);
    }

    void cross(const Vec3d& a, const Vec3d& b) {
      x = a.y * b.z - a.z * b.y;
      y = a.z * b.x - a.x * b.z;
      z = a.x * b.y - a.y * b.x;
    }

    double lengthSqr() const {
      return x * x + y * y + z * z;
    }

    double length() const {
      return sqrt(x * x + y * y + z * z);
    }

    void scale(double s) {
      x *= s;
      y *= s;
      z *= s;
    }

    void scaled(const Vec3d& v, double s) {
      x = v.x * s;
      y = v.y * s;
      z = v.z * s;
    }

    void addScaled(const Vec3d& v, double s) {
      x += v.x * s;
      y += v.y * s;
      z += v.z * s;
    }

    void lineCombin(const Vec3d& a, const Vec3d& b, double t) {
      x = a.x + b.x * t;
      y = a.y + b.y * t;
      z = a.z + b.z * t;
    }

    void lineCombin2(const Vec3d& a, double ta, const Vec3d& b, double tb) {
      x = a.x * ta + b.x * tb;
      y = a.y * ta + b.y * tb;
      z = a.z * ta + b.z * tb;
    }
  };
}


bool equal(const hess::Vec3d& a, const hess::Vec3d& b);
bool is_zero(const hess::Vec3d& v);
void random_vector(hess::Vec3d& a);

