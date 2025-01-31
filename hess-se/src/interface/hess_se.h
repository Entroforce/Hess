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

#ifdef __cplusplus
extern "C" {
#endif
    
  void* hessSECreate();
  
  int hessSEAddMol(void *se, void *atoms, int optMode, int *indices);
  
  int hessSERun(void *se);
  
  int hessSEGetEnergy(void *se,double *energy);
  
  int hessSEGetMol(void *se, int idx, void ** hessMol);
  
#ifdef __cplusplus
}
#endif

  
  
  