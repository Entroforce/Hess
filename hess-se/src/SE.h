/**************************************************************************
 * This file is part of the Hess project
 * Copyright (C) 2024 Entroforce LLC
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

#ifndef SE_H
#define SE_H

#include "model/molecule.h"

struct MopacMol {
    hess::Molecule * atoms;
    int *molMap;
    double *gradient;
    double energy;
    int optMode;
};

class SE {
public:
    SE();
    int addMolecule(hess::Molecule * atoms, int optMode, int *idx);
    virtual ~SE();
    int run();
    int getEnergy(double *energy);
    int getMol(int idx, hess::Molecule ** mol);
    std::string getErrorMessage();
    
    
private:
    MopacMol **mols;
    int numMols;
    std::string error_message;
    double energy;
};

#endif /* SE_H */

