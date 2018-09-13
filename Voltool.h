/***************************************************************************
 *cr
 *cr            (C) Copyright 2007-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: Voltool.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2018/09/12 15:14:14 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  General volumetric data processing routines, particularly supporting MDFF
 *
 ***************************************************************************/

#ifndef VOLTOOL_H
#define VOLTOOL_H

#include <stdio.h>
#include "VolumetricData.h"
#include <stdint.h>
#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

static inline int myisnan(float f)
{
      union { float f; uint32_t x; } u = { f };
          return (u.x << 1) > 0xff000000u;
}

// Cubic interpolation used by supersample
inline float cubic_interp(float y0, float y1, float y2, float y3, float mu) {
  
  float mu2 = mu*mu;
  float a0 = y3 - y2 - y0 + y1;
  float a1 = y0 - y1 - a0;
  float a2 = y2 - y0;
  float a3 = y1;

  return (a0*mu*mu2+a1*mu2+a2*mu+a3);
}

//get the cartesian coordinate of a voxel given its x,y,z indices
inline void voxel_coord(int x, int y, int z, float &gx, float &gy, float &gz, VolumetricData *vol){
  float xdelta[3], ydelta[3], zdelta[3];
  vol->cell_axes(xdelta, ydelta, zdelta);
  
  gx = vol->origin[0] + (x * xdelta[0]) + (y * ydelta[0]) + (z * zdelta[0]);
  gy = vol->origin[1] + (x * xdelta[1]) + (y * ydelta[1]) + (z * zdelta[1]);
  gz = vol->origin[2] + (x * xdelta[2]) + (y * ydelta[2]) + (z * zdelta[2]);
 
}

//get the cartesian coordinate of a voxel given its 1D index 
inline void voxel_coord(int i, float &x, float &y, float &z, VolumetricData *vol){
  float xdelta[3], ydelta[3], zdelta[3];
  vol->cell_axes(xdelta, ydelta, zdelta);
  int xsize = vol->xsize;
  int ysize = vol->ysize;
  
  int gz = i / (ysize*xsize);
  int gy = (i / xsize) % ysize;
  int gx = i % xsize;

  x = vol->origin[0] + (gx * xdelta[0]) + (gy * ydelta[0]) + (gz * zdelta[0]);
  y = vol->origin[1] + (gx * xdelta[1]) + (gy * ydelta[1]) + (gz * zdelta[1]);
  z = vol->origin[2] + (gx * xdelta[2]) + (gy * ydelta[2]) + (gz * zdelta[2]);
 
}

//--unary ops--

//adds or removes voxels in the given axis directions
void pad(VolumetricData *vol, int padxm, int padxp, int padym, int padyp, int padzm, int padzp);
//crops a volumetric data to a given set of cartesian coordinates
void crop(VolumetricData *vol, double crop_minx, double crop_miny, double crop_minz, double crop_maxx, double crop_maxy, double crop_maxz);
//clamp out of range voxel values
void clamp(VolumetricData *vol, float min_value, float max_value);
//scales voxel data by given amount
void scale_by(VolumetricData *vol, float ff);
//adds to voxel data a given amount
void scalar_add(VolumetricData *vol, float ff);
//fits voxel data to a given range
void fit_to_range(VolumetricData *vol, float min_value, float max_value);
//dowmnsample voxels by 2 (x8 total reduction)
void downsample(VolumetricData *vol);
//supersample voxels by 2 (x8 total increase)
void supersample(VolumetricData *vol);
// Transform map to a sigma scale, so that isovalues in VMD correspond
// to number of sigmas above the mean
void sigma_scale(VolumetricData *vol);
// Makes a mask out of a map
// i.e. all values > 0 are set to 1
void binmask(VolumetricData *vol);
//Gaussian blurring (as a 3D convolution)
void gauss3d(VolumetricData *vol, double sigma);
// Fast Gaussian blur that takes advantage of the fact that the dimensions are separable.
void gauss1d(VolumetricData *vol, double sigma);


//--binary ops--
//adds voxel values from two datasets together
void add(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
//subtracts voxel values of mapB from mapA
void subtract(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
//multiplies voxel values of two datasets together
void multiply(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
//creates a new volumetricdata whose values are the average of two given maps at each location
void average(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);

//--volumetric data utilities--
//calculate the centor of mass of a volumetricdata
void vol_com(VolumetricData *vol, float *com);
//move the center of mass of a given volumetricdata to a specific cartesian coordinate
void vol_moveto(VolumetricData *vol, float *com, float *pos);
//apply a matrix transformation to a volumetricdata
void vol_move(VolumetricData *vol,  float *mat);
//create a new volumetricdata from the union of two datasets
void init_from_union(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol);
//create a new volumetricdata from the intersection of two datasets
void init_from_intersection(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol);

#endif

