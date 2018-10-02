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

#include <vector>
#include <utility>
#include <iostream>
#include <sstream>
#include <fstream>
#include "tinyply.h"
#include "ply_io.h"

RequestedProperties::RequestedProperties(std::shared_ptr<tinyply::PlyData> *_tinyPlyDataPtr, std::string _elementName, std::initializer_list<std::string>__properties) {
  tinyPlyDataPtr=_tinyPlyDataPtr;
  elementName=_elementName;
  _properties=__properties;
  for (auto prop : _properties) {
    properties.push_back(std::pair<std::string,bool>(prop,false));
  }
}

int ply_open(
  const char *filename,
  tinyply::PlyFile &file, // output file
  std::vector<RequestedProperties> &requestedPropertiesList,
  bool requestOtherProperties,   // needed if output file is the input file
  bool verbose=false
) {

  int total=0;
  std::ifstream ss(filename, std::ios::binary);
  if (verbose) std::cerr << "Opening " << filename << " ..." << std::endl;
  if (ss.fail()) throw std::runtime_error(std::string("failed to open ") + filename);

  file.parse_header(ss);

  for (auto e : file.get_elements()) {
    if (verbose) std::cerr << "element - " << e.name << " (" << e.size << ")" << std::endl;
    for (auto p : e.properties) {
      bool found=false;
      if (verbose) std::cerr << "\tproperty - " << p.name << " (" << tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
      for(auto req : requestedPropertiesList) {
        if (req.tinyPlyDataPtr[0]->count) continue;
        if (e.name==std::string(req.elementName)) {
          size_t count=0;
          for (std::pair<std::string, bool> prop: req.properties) {
            if (p.name==std::string(prop.first)) {
              prop.second=true;
              found=true;
            }
            if (prop.second) ++count;
          }
          if (count==req.properties.size()) {
            if (true /*verbose*/) {
              std::cerr << "extracted ";
              int i=req._properties.size();
              for (auto name : req._properties) {
                std::cerr << name << ((--i)?", ":"");
              }
            }
            *req.tinyPlyDataPtr=file.request_properties_from_element(e.name.c_str(), req._properties);
            ++total;
          }
          if (found) break;
        }
      }
      if (found) continue;
      (void)file.request_properties_from_element(e.name.c_str(), { p.name.c_str() });
    }
  }

  return total;

}
