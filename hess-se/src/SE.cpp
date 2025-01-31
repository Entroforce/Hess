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

#include "SE.h"
#include "scoring_function/constants.h"
#include <string>
#include <vector>

extern "C" double run_mopac(double *coord, int *labels, char **txtatm, int natoms, const char *keywrd,
                            double *geo_int, int *loc, int *lopt, int *na, int *nb,
                            int *nc, int nvar, double *xparam, double *coords, double *q, double *escf, double *gradient);


SE::SE() {
    this->numMols = 0;
    this->mols = nullptr;
    this->error_message = "";
    this->energy = 0;
}

int SE::addMolecule(hess::Molecule * atoms, int optMode, int *idx){   
    if(atoms == nullptr) {
        this->error_message = "Empty molecule object passed!";
        return -1;
    }
    
    if(atoms->get_atoms_count() < 1) {
        this->error_message = "No atoms in a passed molecule!";
        return -1;
    }
    
    this->numMols++;
    if(this->mols == nullptr) {
        this->mols = (MopacMol **) malloc(sizeof(MopacMol *));
    } else {
        this->mols = (MopacMol **) realloc(this->mols, this->numMols * sizeof(MopacMol *));
    }
    
    MopacMol * mpMol = new MopacMol();
    mpMol->optMode = optMode;
    mpMol->atoms = atoms;
    mpMol->molMap = (int *) malloc(atoms->get_atoms_count() * sizeof(int));
    
    int i = 0;
    
    for (int it = atoms->vertexBegin(); it != atoms->vertexEnd(); it = atoms->vertexNext(it)) {
        mpMol->molMap[i] = it;
        i++;
    }
    
    this->mols[this->numMols - 1] = mpMol;
    
    *idx = this->numMols - 1;

    return 0;
}


SE::~SE() {
    for(int i = 0; i < numMols; i++) {
        MopacMol *mop = this->mols[i];
        free(mop->atoms);
        free(mop->gradient);
        free(mop->molMap);
        free(mop);
    }
    free(this->mols);
    return;
}

int SE::getEnergy(double *energy){
    *energy = this->energy;
    return 0;
}

int SE::run() {
    int numberOfAtoms = 0;
    std::cout << "hrtt" <<std::endl;
    
    for(int i = 0; i < this->numMols; i++){
        if(this->mols[i] == nullptr)
            continue;
        
        MopacMol *mop = this->mols[i];
        if(mop->atoms == nullptr) {
            this->error_message = "Molecule is not specified!";
            return -1;
        }
        
        numberOfAtoms += mop->atoms->get_atoms_count();
    }
   
    vector<double> geo(3 * numberOfAtoms, 0);
    vector<double> coord(3 * numberOfAtoms, 0);
    vector<int> labels(numberOfAtoms, 0);
    vector<int> lopt(3 * numberOfAtoms, 0);
    vector<int> na(numberOfAtoms, 0);
    vector<int> nb(numberOfAtoms, 0);
    vector<int> nc(numberOfAtoms, 0);
    vector<char *> txtatm(numberOfAtoms);
    
    int idx = 0;
    int nvar = 0;
    int isSCF = 1;
    for(int i = 0; i < this->numMols; i++){
        
        int kk = 0;
        MopacMol *mop = this->mols[i];
        
        if(mop->optMode != 0)
            isSCF = 0;
        for (int it = mop->atoms->vertexBegin(); it != mop->atoms->vertexEnd(); it = mop->atoms->vertexNext(it)) {
            labels[idx] = mop->atoms->get_atom(it)->num;
            geo[3 * idx] = mop->atoms->get_atom(it)->x;
            geo[3 * idx + 1] = mop->atoms->get_atom(it)->y;
            geo[3 * idx + 2] = mop->atoms->get_atom(it)->z;
                        
            coord[3 * idx] = mop->atoms->get_atom(it)->x;
            coord[3 * idx + 1] = mop->atoms->get_atom(it)->y;
            coord[3 * idx + 2] = mop->atoms->get_atom(it)->z;
            lopt[3 * idx] = (mop->optMode == 1 || (mop->optMode == 2 && (labels[idx] == 1))) ? 1 : 0;
            lopt[3 * idx + 1] = (mop->optMode == 1 || (mop->optMode == 2 && (labels[idx] == 1))) ? 1 : 0;
            lopt[3 * idx + 2] = (mop->optMode == 1 || (mop->optMode == 2 && (labels[idx] == 1))) ? 1 : 0;
            na[idx] = 0;
            nb[idx] = 0;
            nc[idx] = 0;
            txtatm[idx] = this->mols[i]->atoms->txtatm[kk];
            kk++;
            idx++;
        }
    }
    
    for(int i = 0; i < 3 * numberOfAtoms; i++){
        if(lopt[i] != 0){
            nvar++;
        }
    }
    
    vector<int> loc(2 * nvar, 0);
    vector<double> xparam(nvar, 0);

    if(nvar > 0) {
        int k = 0;
        for(int i = 0; i < numberOfAtoms; i++){
            for(int j = 0; j < 3; j++){
                if(lopt[3 * i + j] == 1){
                    loc[2 * k] = i + 1;
                    loc[2 * k + 1] = j + 1;
                    xparam[k] = geo[3 * i + j];
                    k += 1;
                }
            }
        }
    }
    
    vector<double> coords(3 * numberOfAtoms, 0);
    vector<double> q(numberOfAtoms, 0);
        
    std::string get_charge_keyword = " PM7 MOZYME CHARGES GEO-OK";
    double computed_charge = run_mopac(&coord[0], &labels[0], &txtatm[0], numberOfAtoms, get_charge_keyword.c_str(), &geo[0], &loc[0], &lopt[0], &na[0], &nb[0], &nc[0], nvar, &xparam[0], &coords[0], &q[0], nullptr, nullptr);
    
    std::string keywords;
    if(isSCF == 0)
        keywords = " PM7 PL LBFGS MOZYME EPS=78.4 CHARGE=" + std::to_string(int(computed_charge)) + " SINGLET ";
    else 
        keywords = " PM7 1SCF MOZYME EPS=78.4 CHARGE=" + std::to_string(int(computed_charge)) + " SINGLET ";

    double escf = 0;

    this->energy = run_mopac(&coord[0], &labels[0], &txtatm[0], numberOfAtoms, keywords.c_str(), &geo[0], &loc[0], &lopt[0], &na[0], &nb[0], &nc[0], nvar, &xparam[0], &coords[0], &q[0], &escf, nullptr);

    int ii = 0;
    for(int i = 0; i < this->numMols; i++){
        MopacMol *mop = this->mols[i];
        int numAtoms = mop->atoms->get_atoms_count();
        for(int idx = 0; idx < numAtoms; idx++){
            int it = mop->molMap[idx];
            mop->atoms->setAtomXyz(it, coords[3 * ii], coords[3 * ii + 1], coords[3 * ii + 2]);
            ii++;
        }
    }
     
    return 0;
}

int SE::getMol(int idx, hess::Molecule ** mol) {
    if(idx < 0 || idx >= this->numMols) {
        this->error_message = "Input molecule index out of range!";
        return -1;
    }
    
    if(this->mols[idx] == nullptr) {
        this->error_message = "Molecule with the input index doesn't exist!";
        return -1;
    }
    
    *mol = (this->mols[idx])->atoms;
    
    return 0;
}

std::string SE::getErrorMessage() {
    return this->error_message;
}
