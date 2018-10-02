/*
* Copyright (c) 2018 ALSENET SA
*
* Author(s):
*
*      Luc Deschenaux <luc.deschenaux@freesurf.ch>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#ifndef __PLY_IO__
#define __PLY_IO__

#include <vector>
#include <utility>
#include "tinyply.h"

class RequestedProperties {
public:
  std::string elementName;
  std::shared_ptr<tinyply::PlyData> *tinyPlyDataPtr;
  std::initializer_list<std::string>_properties;
  std::vector<std::pair<std::string, bool> > properties;
  RequestedProperties(std::shared_ptr<tinyply::PlyData> *tinyPlyDataPtr, std::string elementName, std::initializer_list<std::string>_properties);
};


/**
  open a ply file using tinyply and request properties from elements
  @param filename the ply file name
  @param requestedProperties what to extract and where to store the pointers
  @param requestOtherProperties needed if output file is the input file
  @param verbose
  @return status
*/
int ply_open(
  const char *filename,
  tinyply::PlyFile &file, // output file
  std::vector<RequestedProperties> requestedPropertiesList,
  bool requestOtherProperties,   // needed if output file is the input file
  bool verbose=false
);

#endif
