# regply

# LICENSE

Licensed under AGPLv3 or later, see Copyright section below.

# Description

Align a ply file using the point-pairs specified in two other ply files.

# Dependencies

You need at least:

* libpcl-1.8
* libeigen3
* libboost-system

You can install the depencies for example with:
```
sudo apt-get install libpcl-dev libboost-dev libproj-dev
```

## Build and Install

If CloudCompare is not yet installed on your system, build and install it, for example running the following commands in the base directory:
```
cd cloudcompare
mkdir build
cmake .. -DCMAKE_BUILD_TYPE=Release
make && sudo checkinstall --pkgname cloudcompare
```

Then build and install regply, for example running the following commandes in the base directory:
```
mkdir build
cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make && sudo checkinstall --pkgname regply
```

## Synopsis
```
Usage: regply <options>
Options:
  -r|reference <filename>         reference points
  -c|correspondences <filename>   control points to be aligned
  -f|fixed-scale                  do not adjust scale
  -t|transform <filename>         optional: cloud to be transformed using resulting matrix
  -o|output <filename>            optional: output file name
```

## Example

```
regply --reference ref.ply --correspondences control-points.ply --transform cloud.ply --output cloud-aligned.ply
```

## COPYRIGHT

```
 Copyright (c) 2018 ALSENET SA

 Author(s):

      Luc Deschenaux <luc.deschenaux@freesurf.ch>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
```

