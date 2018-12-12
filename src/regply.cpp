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


#include <Eigen/Geometry>

#include <iostream>
#include <proj_api.h>
#include <algorithm>
#include <math.h>
#include <typeinfo>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <getopt.h>

#include "PointCloud.h"
#include "PointProjectionTools.h"
#include "ScalarFieldTools.h"
#include "RegistrationTools.h"
#include "tinyply.h"
#include "ply_io.h"
#include "regply.h"

using namespace CCLib;

struct float3 { float x, y, z; };
struct double3 { double x, y, z; };

int registration(char *ref_filename, char *cor_filename, char *align_filename, char *out_filename, bool fixedScale, bool shift);

static int fixedScale;
static int shiftCoordinates;
static int verbose;
char *ref_filename;
char *cor_filename;
char *align_filename;
char *out_filename;
char *matrix;
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
  std::cerr << "  -s|shift                        shift coordinates" << std::endl;
  std::cerr << "  -a|align <filename>             optional: cloud to be aligned using resulting matrix" << std::endl;
  std::cerr << "  -o|output <filename>            optional: output file name" << std::endl;
  std::cerr << "  -v|verbose" << std::endl;
  exit(1);
}


int main(int argc, char **argv) {
  appName=argv[0];
  int c;
  while(1) {
    static struct option long_options[] = {
      {"fixed-scale", no_argument, &fixedScale, 1},
      {"shift", no_argument, &shiftCoordinates, 1},
      {"verbose", no_argument, &verbose, 1},
      {"reference", required_argument, 0, 'r'},
      {"correspondences", required_argument, 0, 'c'},
      {"align", required_argument, 0, 'a'},
      {"output", required_argument, 0, 'o'},
      {"help", required_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

    int option_index = 0;

    c = getopt_long (argc, argv, "fr:c:a:o:hvs", long_options, &option_index);

    if (c == -1)
      break;

    switch(c){
      case 'r':
        ref_filename=optarg;
        break;

      case 'c':
        cor_filename=optarg;
        break;

      case 'a':
        align_filename=optarg;
        break;

      case 'o':
        out_filename=optarg;
        break;

      case 'f':
        fixedScale=1;
        break;

      case 's':
        shiftCoordinates=1;
        break;

      case 'v':
        verbose=1;
        break;

      default:
        usage();
        break;
    }
  }

  if (!ref_filename) {
    std::cerr << "You must specify the reference points filename with -r" << std::endl;
    usage();
  }

  if (!cor_filename) {
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

  try {
    return registration(ref_filename,cor_filename,align_filename,out_filename,fixedScale,shiftCoordinates);
  } catch(const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
}

template<class T>
void getIdealShift(T *points, size_t count, T *shift) {
  if (!count) {
    shift->x=shift->y=shift->z=0;
    return;
  }

  typeof(points->x) minx,miny,minz;
  minx=std::numeric_limits<typeof(points->x)>::max();
  miny=std::numeric_limits<typeof(points->x)>::max();
  minz=std::numeric_limits<typeof(points->x)>::max();

  for (size_t i = 0; i < count; ++i) {
    T *point=points+i;
    minx=std::min(minx,point->x);
    miny=std::min(miny,point->y);
    minz=std::min(minz,point->z);
   }

   shift->x=static_cast<long int>(minx/100)*100.0;
   shift->y=static_cast<long int>(miny/100)*100.0;
   shift->z=static_cast<long int>(minz/100)*100.0;
}

template<class T, class U>
void fillCloud(PointCloud &P, PointCloud &X, T cor, U ref, size_t count) {
   for (size_t i = 0; i < count; ++i) {
     P.addPoint(CCVector3(cor[i].x, cor[i].y, cor[i].z));
     X.addPoint(CCVector3(ref[i].x, ref[i].y, ref[i].z));
   }
}

template<class T, class U>
void fillCloud(PointCloud &P, PointCloud &X, T cor, U ref, size_t count, T cor_offset, U ref_offset) {
   for (size_t i = 0; i < count; ++i) {
     P.addPoint(CCVector3(cor[i].x-cor_offset->x, cor[i].y-cor_offset->y, cor[i].z-cor_offset->z));
     X.addPoint(CCVector3(ref[i].x-ref_offset->x, ref[i].y-ref_offset->y, ref[i].z-ref_offset->z));
   }
}

void streamTrans(std::ostream &stream, PointProjectionTools::Transformation &trans, double rms) {
  // print the transformation matrix
  stream << "# Final RMS: " << rms << std::endl;
  if (!fixedScale) stream << "# Scale: " << trans.s  << " (already integrated in matrix below) " << std::endl;
  stream << std::fixed;
  stream << trans.R.getValue(0,0)*trans.s << ", " << trans.R.getValue(0,1)*trans.s << ", " << trans.R.getValue(0,2)*trans.s << ", " << trans.T.x << std::endl;
  stream << trans.R.getValue(1,0)*trans.s << ", " << trans.R.getValue(1,1)*trans.s << ", " << trans.R.getValue(1,2)*trans.s << ", " << trans.T.y << std::endl;
  stream << trans.R.getValue(2,0)*trans.s << ", " << trans.R.getValue(2,1)*trans.s << ", " << trans.R.getValue(2,2)*trans.s << ", " << trans.T.z << std::endl;
  stream << 0.0 << ", " << 0.0 << ", " << 0.0 << ", " << 1.0 << std::endl;
}

template<class T, class U, class V>
void apply_matrix(T points, size_t count, Eigen::Affine3d eigen_trans, U cor_offset, V ref_offset) {
  Eigen::Vector4d v(0,0,0,1);

  double cox=cor_offset->x;
  double coy=cor_offset->y;
  double coz=cor_offset->z;
  double rox=ref_offset->x;
  double roy=ref_offset->y;
  double roz=ref_offset->z;

  if (typeid(points->x)==typeid(float)) {
    for (size_t i=0; i<count; ++i) {
      T p=points+i;
      v(0)=static_cast<double>(p->x)-cox;
      v(1)=static_cast<double>(p->y)-coy;
      v(2)=static_cast<double>(p->z)-coz;
      Eigen::Vector4d w=eigen_trans*v;
      p->x=static_cast<float>(w(0)+rox);
      p->y=static_cast<float>(w(1)+roy);
      p->z=static_cast<float>(w(2)+roz);
    }
  } else {
    for (size_t i=0; i<count; ++i) {
      T p=points+i;
      v(0)=p->x-cox;
      v(1)=p->y-coy;
      v(2)=p->z-coz;
      Eigen::Vector4d w=eigen_trans*v;
      p->x=w(0)+rox;
      p->y=w(1)+roy;
      p->z=w(2)+roz;
    }
  }
}

template<class T>
void apply_matrix(T points, size_t count, Eigen::Affine3d eigen_trans) {
  Eigen::Vector4d v(0,0,0,1);
  if (typeid(points->x)==typeid(float)) {
    for (size_t i=0; i<count; ++i) {
      T p=points+i;
      v(0)=static_cast<double>(p->x);
      v(1)=static_cast<double>(p->y);
      v(2)=static_cast<double>(p->z);
      Eigen::Vector4d w=eigen_trans*v;
      p->x=static_cast<float>(w(0));
      p->y=static_cast<float>(w(1));
      p->z=static_cast<float>(w(2));
    }
  } else {
    for (size_t i=0; i<count; ++i) {
      T p=points+i;
      v(0)=p->x;
      v(1)=p->y;
      v(2)=p->z;
      Eigen::Vector4d w=eigen_trans*v;
      p->x=w(0);
      p->y=w(1);
      p->z=w(2);
    }
  }
}

int registration(char *ref_filename, char *cor_filename, char *align_filename, char *out_filename, bool fixedScale, bool shiftCoordinates) {

  tinyply::PlyFile ref_file;
  tinyply::PlyFile cor_file;

  std::shared_ptr<tinyply::PlyData> ref_vertices=0;
  std::vector<RequestedProperties> ref_requestList;
  ref_requestList.push_back(RequestedProperties(&ref_vertices,"vertex",{"x","y","z"}));

  std::shared_ptr<tinyply::PlyData> cor_vertices=0;
  std::vector<RequestedProperties> cor_requestList;
  cor_requestList.push_back(RequestedProperties(&cor_vertices,"vertex",{"x","y","z"}));

  ply_read(ref_filename,ref_file,ref_requestList,verbose);
  ply_read(cor_filename,cor_file,cor_requestList,verbose);

  if (!ref_vertices) throw std::runtime_error(std::string("failed to extract vertice properties from ") + ref_filename);
  if (!cor_vertices) throw std::runtime_error(std::string("failed to extract vertice properties from ") + cor_filename);

  if (ref_vertices->count!=cor_vertices->count) {
    throw std::runtime_error(std::string("number of points must be equal in both files"));
  }

  size_t count=ref_vertices->count;

  // fill CC pointclouds
  PointCloud P,X;
  void *ref_points_buf=ref_vertices->buffer.get();
  void *cor_points_buf=cor_vertices->buffer.get();

  float3 ref_shift_f;
  float3 cor_shift_f;
  double3 ref_shift_d;
  double3 cor_shift_d;

  if (shiftCoordinates) {
    if (ref_vertices->t==tinyply::Type::FLOAT32) {
      getIdealShift((float3*)ref_points_buf, count, &ref_shift_f);
      if (cor_vertices->t==tinyply::Type::FLOAT32) {
        getIdealShift((float3*)cor_points_buf, count, &cor_shift_f);
        fillCloud(P, X, (float3*)cor_points_buf, (float3*)ref_points_buf, count, &cor_shift_f, &ref_shift_f);
      } else {
        getIdealShift((double3*)cor_points_buf, count, &cor_shift_d);
        fillCloud(P, X, (double3*)cor_points_buf, (float3*)ref_points_buf, count, &cor_shift_d, &ref_shift_f);
      }
    } else {
      getIdealShift((double3*)ref_points_buf, count, &ref_shift_d);
      if (cor_vertices->t==tinyply::Type::FLOAT32) {
        getIdealShift((float3*)cor_points_buf, count, &cor_shift_f);
        fillCloud(P, X, (float3*)cor_points_buf, (double3*)ref_points_buf, count, &cor_shift_f, &ref_shift_d);
      } else {
        getIdealShift((double3*)cor_points_buf, count, &cor_shift_d);
        fillCloud(P, X, (double3*)cor_points_buf, (double3*)ref_points_buf, count, &cor_shift_d, &ref_shift_d);
      }
    }
  } else {
    if (ref_vertices->t==tinyply::Type::FLOAT32) {
      if (cor_vertices->t==tinyply::Type::FLOAT32) {
        fillCloud(P, X, (float3*)cor_points_buf, (float3*)ref_points_buf, count);
      } else {
        fillCloud(P, X, (double3*)cor_points_buf, (float3*)ref_points_buf, count);
      }
    } else {
      if (cor_vertices->t==tinyply::Type::FLOAT32) {
        fillCloud(P, X, (float3*)cor_points_buf, (double3*)ref_points_buf, count);
      } else {
        fillCloud(P, X, (double3*)cor_points_buf, (double3*)ref_points_buf, count);
      }
    }
  }

  // CC registration
  PointProjectionTools::Transformation trans;
  double rms;
  if (HornRegistrationTools::FindAbsoluteOrientation((GenericCloud*)&P,(GenericCloud*)&X,trans,fixedScale)) {
    rms=HornRegistrationTools::ComputeRMS((GenericCloud*)&P,(GenericCloud*)&X,trans);
  } else {
    std::cerr << "Registration failed !" << std::endl;
    exit (1);
  }

  if (align_filename) {
    // align specified ply
    Eigen::Affine3d eigen_trans;
    eigen_trans(0,0)=static_cast<double>(trans.R.getValue(0,0))*static_cast<double>(trans.s);
    eigen_trans(0,0)=static_cast<double>(trans.R.getValue(0,1))*static_cast<double>(trans.s);
    eigen_trans(0,0)=static_cast<double>(trans.R.getValue(0,2))*static_cast<double>(trans.s);
    eigen_trans(0,3)=static_cast<double>(trans.T.x);
    eigen_trans(1,0)=static_cast<double>(trans.R.getValue(1,0))*static_cast<double>(trans.s);
    eigen_trans(1,0)=static_cast<double>(trans.R.getValue(1,1))*static_cast<double>(trans.s);
    eigen_trans(1,0)=static_cast<double>(trans.R.getValue(1,2))*static_cast<double>(trans.s);
    eigen_trans(1,3)=static_cast<double>(trans.T.y);
    eigen_trans(2,0)=static_cast<double>(trans.R.getValue(2,0))*static_cast<double>(trans.s);
    eigen_trans(2,0)=static_cast<double>(trans.R.getValue(2,1))*static_cast<double>(trans.s);
    eigen_trans(2,0)=static_cast<double>(trans.R.getValue(2,2))*static_cast<double>(trans.s);
    eigen_trans(2,3)=static_cast<double>(trans.T.z);

    std::shared_ptr<tinyply::PlyData> align_vertices=0;
    std::vector<RequestedProperties> align_requestList;
    align_requestList.push_back(RequestedProperties(&align_vertices,"vertex",{"x","y","z"}));

    tinyply::PlyFile align_file;
    ply_read(align_filename,align_file,align_requestList,verbose);
    if (!align_vertices) throw std::runtime_error(std::string("failed to extract vertice properties from ") + align_filename);

    std::vector<std::string> &comments=align_file.get_comments();
    comments.push_back(std::string("regply final RMS: ")+std::to_string(rms));

    void *points_buf=align_vertices->buffer.get();
    Eigen::Vector4d v(0,0,0,1);

    if (shiftCoordinates) {
      if (cor_vertices->t==tinyply::Type::FLOAT32) {
        cor_shift_d.x=static_cast<double>(cor_shift_f.x);
        cor_shift_d.y=static_cast<double>(cor_shift_f.y);
        cor_shift_d.z=static_cast<double>(cor_shift_f.z);
      }
      if (ref_vertices->t==tinyply::Type::FLOAT32) {
        ref_shift_d.x=static_cast<double>(ref_shift_f.x);
        ref_shift_d.y=static_cast<double>(ref_shift_f.y);
        ref_shift_d.z=static_cast<double>(ref_shift_f.z);
      }

      if (align_vertices->t==tinyply::Type::FLOAT32) {
        apply_matrix(
          (float3*)points_buf,
          align_vertices->count,
          eigen_trans,
          &cor_shift_d,
          &ref_shift_d
        );
      } else {
        apply_matrix(
          (double3*)points_buf,
          align_vertices->count,
          eigen_trans,
          &cor_shift_d,
          &ref_shift_d
        );
      }
    } else {
      if (align_vertices->t==tinyply::Type::FLOAT32) {
        apply_matrix(
          (float3*)points_buf,
          align_vertices->count,
          eigen_trans
        );
      } else {
        apply_matrix(
          (double3*)points_buf,
          align_vertices->count,
          eigen_trans
        );
      }
    }

    if (out_filename) {
      // save to specified output file
      std::filebuf fb_ascii;
      fb_ascii.open(out_filename, std::ios::out);
      std::ostream outstream_ascii(&fb_ascii);
      if (outstream_ascii.fail()) throw std::runtime_error(std::string("failed to open ") + out_filename + " for writing");
      outstream_ascii << std::fixed;
      align_file.write(outstream_ascii,false);

    } else {
      // or print to stdout
      std::cout << std::fixed;
      align_file.write(std::cout,false);
    }

  }

  if (!align_filename) {
    if (out_filename) {
      // save to specified output file
      std::filebuf fb_ascii;
      fb_ascii.open(out_filename, std::ios::out);
      std::ostream outstream_ascii(&fb_ascii);
      if (outstream_ascii.fail()) throw std::runtime_error(std::string("failed to open ") + out_filename + " for writing");
      streamTrans(outstream_ascii,trans,rms);

    } else {
      // or print to stdout
      streamTrans(std::cout,trans,rms);
    }

  }

  return 0;

}
