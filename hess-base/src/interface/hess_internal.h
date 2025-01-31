/**************************************************************************
 * This file is part of the Hess project
 * Copyright (C) 2023-2025 Entroforce LLC
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
#include "parser/parser.h"
#include "empirical-docking/optimizable_molecule.h"

struct HessObject {

  HessObject(hess::Molecule* _mol = nullptr, Optimizable_molecule* _opt_mol = nullptr, hess::Parser* _parser = nullptr) : mol(_mol), opt_mol(_opt_mol), parser(_parser) {
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

  hess::Molecule* mol;
  Optimizable_molecule* opt_mol;
  hess::Parser* parser;
};

extern thread_local std::string hess_error_message;