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

#include "PointCloud.h"
#include "PointProjectionTools.h"
#include "ScalarFieldTools.h"
#include "RegistrationTools.h"
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/common/transformation_from_correspondences.h>
#include <pcl/common/transforms.h>

#include <Eigen/Geometry>

#include <iostream>
#include <proj_api.h>
#include <algorithm>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "regply.h"

using namespace CCLib;

int registration(char *ref_filename, char *cor_filename, char *align_filename, char *out_filename, bool fixedScale);

static int fixedScale=0;
char *reference;
char *correspondences;
char *transform;
char *output;
char *appName;

void version() {
  std::cerr << appName << " " \
    << " Version " << regply_VERSION_MAJOR << "." << regply_VERSION_MINOR << "." << regply_VERSION_PATCH \
    << ", branch " << regply_GIT_BRANCH << ", commit " << regply_GIT_COMMIT << std::endl;
}

void usage() {
  version();
  std::cerr << "Usage: " << appName << " <options>" << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << "  -r|reference <filename>         reference points" << std::endl;
  std::cerr << "  -c|correspondences <filename>   control points to be aligned" << std::endl;
  std::cerr << "  -f|fixed-scale                  do not adjust scale" << std::endl;
  std::cerr << "  -t|transform <filename>         optional: cloud to be transformed using resulting matrix" << std::endl;
  std::cerr << "  -o|output <filename>            optional: output file name" << std::endl;
  exit(1);
}


int main(int argc, char **argv) {
  appName=argv[0];
  int c;
  while(1) {
    static struct option long_options[] = {
      {"fixed-scale", no_argument, &fixedScale, 1},
      {"reference", required_argument, 0, 'r'},
      {"correspondences", required_argument, 0, 'c'},
      {"transform", required_argument, 0, 't'},
      {"output", required_argument, 0, 'o'},
      {"help", required_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "fr:c:t:o:h", long_options, &option_index);

    if (c == -1)
      break;

    switch(c){
      case 'r':
        reference=optarg;
        break;

      case 'c':
        correspondences=optarg;
        break;

      case 't':
        transform=optarg;
        break;

      case 'o':
        output=optarg;
        break;

      case 'f':
        fixedScale=1;
        break;

      default:
        usage();
        break;
    }
  }

  if (!reference) {
    std::cerr << "You must specify the reference points filename with -r" << std::endl;
    usage();
  }

  if (!correspondences) {
    std::cerr << "You must specify the correspondences points filename with -c" << std::endl;
    usage();
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc) {
      std::cerr << "invalid arguments:";
      while (optind < argc) {
        std::cerr << " " << argv[optind++];
      }
      std::cerr << std::endl;
      std::cerr << std::endl;
      usage();
    }

  return registration(reference,correspondences,transform,output,fixedScale);
}

int registration(char *ref_filename, char *cor_filename, char *align_filename, char *out_filename, bool fixedScale) {

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PointCloud<pcl::PointXYZ>::Ptr refCloud (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PLYReader Reader;
  Reader.read(ref_filename, *refCloud);
  Reader.read(cor_filename, *cloud);
  PointCloud P,X;
  PointProjectionTools::Transformation trans;
  Eigen::Affine3f eigen_trans;

  pcl::TransformationFromCorrespondences pcl_transformation_from_correspondences;

  if (cloud->points.size()!=refCloud->points.size()) {
    std::cerr << "number of points must be equal in both files" << std::endl;
    return false;
  }

  for (size_t i = 0; i < cloud->points.size(); ++i) {
    P.addPoint(CCVector3(cloud->points[i].x, cloud->points[i].y, cloud->points[i].z));
    X.addPoint(CCVector3(refCloud->points[i].x, refCloud->points[i].y, refCloud->points[i].z));
  }

  if (HornRegistrationTools::FindAbsoluteOrientation((GenericCloud*)&P,(GenericCloud*)&X,trans,fixedScale)) {
    double rms=HornRegistrationTools::ComputeRMS((GenericCloud*)&P,(GenericCloud*)&X,trans);

    std::cout << "Final RMS: " << rms << std::endl;

    std::cout << "-------------------" << std::endl;
    std::cout << "Transformation matrix" << std::endl;

    std::cout << trans.R.getValue(0,0)*trans.s << " " << trans.R.getValue(0,1)*trans.s << " " << trans.R.getValue(0,2)*trans.s << " " /*<< trans.R.getValue(0,3)*trans.s << " " */ << trans.T.x << std::endl;
    std::cout << trans.R.getValue(1,0)*trans.s << " " << trans.R.getValue(1,1)*trans.s << " " << trans.R.getValue(1,2)*trans.s << " " /* << trans.R.getValue(1,3)*trans.s << " " */ << trans.T.y << std::endl;
    std::cout << trans.R.getValue(2,0)*trans.s << " " << trans.R.getValue(2,1)*trans.s << " " << trans.R.getValue(2,2)*trans.s << " " /* << trans.R.getValue(2,3)*trans.s << " " */ << trans.T.z << std::endl;
    std::cout << "0.0 0.0 0.0 1.0" << std::endl;
    std::cout << "-------------------" << std::endl;
    if (!fixedScale) std::cout << "Scale: " << trans.s  << " (already integrated in above matrix) " << std::endl;
    std::cout << std::endl;
  } else {
    std::cerr << "Registration failed !" << std::endl;
    exit (1);
  }

  eigen_trans(0,0)=trans.R.getValue(0,0)*trans.s;
  eigen_trans(0,1)=trans.R.getValue(0,1)*trans.s;
  eigen_trans(0,2)=trans.R.getValue(0,2)*trans.s;
  eigen_trans(0,3)=trans.T.x;
  eigen_trans(1,0)=trans.R.getValue(1,0)*trans.s;
  eigen_trans(1,1)=trans.R.getValue(1,1)*trans.s;
  eigen_trans(1,2)=trans.R.getValue(1,2)*trans.s;
  eigen_trans(1,3)=trans.T.y;
  eigen_trans(2,0)=trans.R.getValue(2,0)*trans.s;
  eigen_trans(2,1)=trans.R.getValue(2,1)*trans.s;
  eigen_trans(2,2)=trans.R.getValue(2,2)*trans.s;
  eigen_trans(2,3)=trans.T.z;
  eigen_trans(3,0)=0;
  eigen_trans(3,1)=0;
  eigen_trans(3,2)=0;
  eigen_trans(3,3)=1;

  if (align_filename) {
    Reader.read(align_filename, *cloud);
  }

  pcl::transformPointCloud(*cloud,*cloud,eigen_trans,false);

  if (out_filename) {
    pcl::io::savePLYFile(out_filename,*cloud,false);

  } else {
    std::cerr << "result:" << std::endl;
    for (size_t i = 0; i < cloud->points.size (); ++i) {
      std::cout << cloud->points[i].x
      << " "    << cloud->points[i].y
      << " "    << cloud->points[i].z
      << std::endl;
    }
  }

  return 0;

}
