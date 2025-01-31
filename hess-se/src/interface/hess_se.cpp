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

#include "interface/hess_internal.h"
#include "SE.h"

extern "C" void * hessSECreate() {
    return (void *) new SE();
}

extern "C" int hessSEAddMol(void *se, void *atoms, int optMode, int *indices) {
    SE* obj = (SE *) se;
    
    HessObject* ho = (HessObject*) atoms;
    if (ho->mol == nullptr) {
      hess_error_message = "The object is not a molecule";
      return -1;
    }
    
    int res = obj->addMolecule(ho->mol, optMode, indices);
    hess_error_message = obj->getErrorMessage();
    
    return res;
}

extern "C" int hessSERun(void *se) {
    SE* obj = (SE *) se;
    
    int res = obj->run();
    hess_error_message = obj->getErrorMessage();
    return res;
}

extern "C" int hessSEGetEnergy(void *se, double *energy) {
    SE* obj = (SE *) se;
    int res = obj->getEnergy(energy);
    hess_error_message = obj->getErrorMessage();
    return res;
}

extern "C" int hessSEGetMol(void *se, int idx, void ** hessMol) {
    SE *obj = (SE *) se;
    hess::Molecule *mol;
    int res = obj->getMol(idx, &mol);
    hess_error_message = obj->getErrorMessage();
    if(res != 0)
        return res;
    *hessMol = (void *)new HessObject(mol, nullptr, nullptr);
    return 0;
}
