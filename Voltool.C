#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h> // FLT_MAX etc
#include "Inform.h"
#include "utilities.h"
//#include "SymbolTable.h"

#include "AtomSel.h"
#include "VMDApp.h"
#include "MoleculeList.h"
#include "VolumetricData.h"
#include "VolMapCreate.h" // volmap_write_dx_file()

#include "CUDAMDFF.h"
#include "MDFF.h"
#include <math.h>
#include <tcl.h>
#include "TclCommands.h"
#include "Measure.h"
#include "MolFilePlugin.h"

#include <iostream>
#include <string>
#include <sstream>
#include "Voltool.h"

/// creates axes, bounding box and allocates data based on 
/// geometrical intersection of A and B
void init_from_intersection(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol) {
  int d;
  
  // Find intersection of A and B
  // The following has been verified for orthog. cells
  // (Does not work for non-orthog cells)
  
  for (d=0; d<3; d++) {
    newvol->origin[d] = MAX(mapA->origin[d], mapB->origin[d]);
    newvol->xaxis[d] = MAX(MIN(mapA->origin[d]+mapA->xaxis[d], mapB->origin[d]+mapB->xaxis[d]), newvol->origin[d]);
    newvol->yaxis[d] = MAX(MIN(mapA->origin[d]+mapA->yaxis[d], mapB->origin[d]+mapB->yaxis[d]), newvol->origin[d]);
    newvol->zaxis[d] = MAX(MIN(mapA->origin[d]+mapA->zaxis[d], mapB->origin[d]+mapB->zaxis[d]), newvol->origin[d]);
  }
    
  vec_sub(newvol->xaxis, newvol->xaxis, newvol->origin);
  vec_sub(newvol->yaxis, newvol->yaxis, newvol->origin);
  vec_sub(newvol->zaxis, newvol->zaxis, newvol->origin);
  
  newvol->xsize = (int) MAX(dot_prod(newvol->xaxis,mapA->xaxis)*mapA->xsize/dot_prod(mapA->xaxis,mapA->xaxis), \
                    dot_prod(newvol->xaxis,mapB->xaxis)*mapB->xsize/dot_prod(mapB->xaxis,mapB->xaxis));
  newvol->ysize = (int) MAX(dot_prod(newvol->yaxis,mapA->yaxis)*mapA->ysize/dot_prod(mapA->yaxis,mapA->yaxis), \
                    dot_prod(newvol->yaxis,mapB->yaxis)*mapB->ysize/dot_prod(mapB->yaxis,mapB->yaxis));
  newvol->zsize = (int) MAX(dot_prod(newvol->zaxis,mapA->zaxis)*mapA->zsize/dot_prod(mapA->zaxis,mapA->zaxis), \
                    dot_prod(newvol->zaxis,mapB->zaxis)*mapB->zsize/dot_prod(mapB->zaxis,mapB->zaxis));
/*   
  for (d=0; d<3; d++) {
    newvol->xdelta[d] = newvol->xaxis[d]/(newvol->xsize-1);
    newvol->ydelta[d] = newvol->yaxis[d]/(newvol->ysize-1);
    newvol->zdelta[d] = newvol->zaxis[d]/(newvol->zsize-1);
  }
 */ 
  // Create map...
  if (newvol->data) delete[] newvol->data;
  newvol->data = new float[newvol->xsize*newvol->ysize*newvol->zsize];
}

/// creates axes, bounding box and allocates data based on 
/// geometrical intersection of A and B
void init_from_intersection(VolumetricData *mapA, VolumetricData *mapB, VolumetricData *newvol) {
  int d;
  
  // Find intersection of A and B
  // The following has been verified for orthog. cells
  // (Does not work for non-orthog cells)
  
  for (d=0; d<3; d++) {
    newvol->origin[d] = MAX(mapA->origin[d], mapB->origin[d]);
    newvol->xaxis[d] = MAX(MIN(mapA->origin[d]+mapA->xaxis[d], mapB->origin[d]+mapB->xaxis[d]), newvol->origin[d]);
    newvol->yaxis[d] = MAX(MIN(mapA->origin[d]+mapA->yaxis[d], mapB->origin[d]+mapB->yaxis[d]), newvol->origin[d]);
    newvol->zaxis[d] = MAX(MIN(mapA->origin[d]+mapA->zaxis[d], mapB->origin[d]+mapB->zaxis[d]), newvol->origin[d]);
  }
    
  vec_sub(newvol->xaxis, newvol->xaxis, newvol->origin);
  vec_sub(newvol->yaxis, newvol->yaxis, newvol->origin);
  vec_sub(newvol->zaxis, newvol->zaxis, newvol->origin);
  
  newvol->xsize = (int) MAX(dot_prod(newvol->xaxis,mapA->xaxis)*mapA->xsize/dot_prod(mapA->xaxis,mapA->xaxis), \
                    dot_prod(newvol->xaxis,mapB->xaxis)*mapB->xsize/dot_prod(mapB->xaxis,mapB->xaxis));
  newvol->ysize = (int) MAX(dot_prod(newvol->yaxis,mapA->yaxis)*mapA->ysize/dot_prod(mapA->yaxis,mapA->yaxis), \
                    dot_prod(newvol->yaxis,mapB->yaxis)*mapB->ysize/dot_prod(mapB->yaxis,mapB->yaxis));
  newvol->zsize = (int) MAX(dot_prod(newvol->zaxis,mapA->zaxis)*mapA->zsize/dot_prod(mapA->zaxis,mapA->zaxis), \
                    dot_prod(newvol->zaxis,mapB->zaxis)*mapB->zsize/dot_prod(mapB->zaxis,mapB->zaxis));
/*   
  for (d=0; d<3; d++) {
    newvol->xdelta[d] = newvol->xaxis[d]/(newvol->xsize-1);
    newvol->ydelta[d] = newvol->yaxis[d]/(newvol->ysize-1);
    newvol->zdelta[d] = newvol->zaxis[d]/(newvol->zsize-1);
  }
 */ 
  // Create map...
  if (newvol->data) delete[] newvol->data;
  newvol->data = new float[newvol->xsize*newvol->ysize*newvol->zsize];
}


/// creates axes, bounding box and allocates data based on 
/// geometrical union of A and B
void init_from_union(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol) {
  // Find union of A and B
  // The following has been verified for orthog. cells
  // (Does not work for non-orthog cells)
  
  vec_zero(newvol->xaxis);
  vec_zero(newvol->yaxis);
  vec_zero(newvol->zaxis);
  
  int d;
  for (d=0; d<3; d++) {
    newvol->origin[d] = MIN(mapA->origin[d], mapB->origin[d]);
  }
  d=0;
  newvol->xaxis[d] = MAX(MAX(mapA->origin[d]+mapA->xaxis[d], mapB->origin[d]+mapB->xaxis[d]), newvol->origin[d]);
  d=1;
  newvol->yaxis[d] = MAX(MAX(mapA->origin[d]+mapA->yaxis[d], mapB->origin[d]+mapB->yaxis[d]), newvol->origin[d]);
  d=2;
  newvol->zaxis[d] = MAX(MAX(mapA->origin[d]+mapA->zaxis[d], mapB->origin[d]+mapB->zaxis[d]), newvol->origin[d]);
  
  newvol->xaxis[0] -= newvol->origin[0];
  newvol->yaxis[1] -= newvol->origin[1];
  newvol->zaxis[2] -= newvol->origin[2];
  
  newvol->xsize = (int) MAX(dot_prod(newvol->xaxis,mapA->xaxis)*mapA->xsize/dot_prod(mapA->xaxis,mapA->xaxis), \
                    dot_prod(newvol->xaxis,mapB->xaxis)*mapB->xsize/dot_prod(mapB->xaxis,mapB->xaxis));
  newvol->ysize = (int) MAX(dot_prod(newvol->yaxis,mapA->yaxis)*mapA->ysize/dot_prod(mapA->yaxis,mapA->yaxis), \
                    dot_prod(newvol->yaxis,mapB->yaxis)*mapB->ysize/dot_prod(mapB->yaxis,mapB->yaxis));
  newvol->zsize = (int) MAX(dot_prod(newvol->zaxis,mapA->zaxis)*mapA->zsize/dot_prod(mapA->zaxis,mapA->zaxis), \
                    dot_prod(newvol->zaxis,mapB->zaxis)*mapB->zsize/dot_prod(mapB->zaxis,mapB->zaxis));
/*  
  for (d=0; d<3; d++) {
    xdelta[d] = xaxis[d]/(xsize-1);
    ydelta[d] = yaxis[d]/(ysize-1);
    zdelta[d] = zaxis[d]/(zsize-1);
  }
  */
  // Create map...
  if (newvol->data) delete[] newvol->data;
  newvol->data = new float[newvol->xsize*newvol->ysize*newvol->zsize];
}



void init_from_identity(VolumetricData *mapA, VolumetricData *newvol) {

  vec_copy(newvol->origin, mapA->origin);
  vec_copy(newvol->xaxis, mapA->xaxis);
  vec_copy(newvol->yaxis, mapA->yaxis);
  vec_copy(newvol->zaxis, mapA->zaxis); 
  
  newvol->xsize = mapA->xsize;
  newvol->ysize = mapA->ysize;
  newvol->zsize = mapA->zsize;
/*    
  int d;
  for (d=0; d<3; d++) {
    xdelta[d] = xaxis[d]/(xsize-1);
    ydelta[d] = yaxis[d]/(ysize-1);
    zdelta[d] = zaxis[d]/(zsize-1);
  }
  */
  // Create map...
  if (newvol->data) delete[] newvol->data;
  newvol->data = new float[newvol->xsize*newvol->ysize*newvol->zsize];
}


/// creates axes, bounding box and allocates data based on 
/// geometrical union of A and B
void init_from_union(VolumetricData *mapA, VolumetricData *mapB, VolumetricData *newvol) {
  // Find union of A and B
  // The following has been verified for orthog. cells
  // (Does not work for non-orthog cells)
  
  vec_zero(newvol->xaxis);
  vec_zero(newvol->yaxis);
  vec_zero(newvol->zaxis);
  
  int d;
  for (d=0; d<3; d++) {
    newvol->origin[d] = MIN(mapA->origin[d], mapB->origin[d]);
  }
  d=0;
  newvol->xaxis[d] = MAX(MAX(mapA->origin[d]+mapA->xaxis[d], mapB->origin[d]+mapB->xaxis[d]), newvol->origin[d]);
  d=1;
  newvol->yaxis[d] = MAX(MAX(mapA->origin[d]+mapA->yaxis[d], mapB->origin[d]+mapB->yaxis[d]), newvol->origin[d]);
  d=2;
  newvol->zaxis[d] = MAX(MAX(mapA->origin[d]+mapA->zaxis[d], mapB->origin[d]+mapB->zaxis[d]), newvol->origin[d]);
  
  newvol->xaxis[0] -= newvol->origin[0];
  newvol->yaxis[1] -= newvol->origin[1];
  newvol->zaxis[2] -= newvol->origin[2];
  
  newvol->xsize = (int) MAX(dot_prod(newvol->xaxis,mapA->xaxis)*mapA->xsize/dot_prod(mapA->xaxis,mapA->xaxis), \
                    dot_prod(newvol->xaxis,mapB->xaxis)*mapB->xsize/dot_prod(mapB->xaxis,mapB->xaxis));
  newvol->ysize = (int) MAX(dot_prod(newvol->yaxis,mapA->yaxis)*mapA->ysize/dot_prod(mapA->yaxis,mapA->yaxis), \
                    dot_prod(newvol->yaxis,mapB->yaxis)*mapB->ysize/dot_prod(mapB->yaxis,mapB->yaxis));
  newvol->zsize = (int) MAX(dot_prod(newvol->zaxis,mapA->zaxis)*mapA->zsize/dot_prod(mapA->zaxis,mapA->zaxis), \
                    dot_prod(newvol->zaxis,mapB->zaxis)*mapB->zsize/dot_prod(mapB->zaxis,mapB->zaxis));
/*  
  for (d=0; d<3; d++) {
    xdelta[d] = xaxis[d]/(xsize-1);
    ydelta[d] = yaxis[d]/(ysize-1);
    zdelta[d] = zaxis[d]/(zsize-1);
  }
  */
  // Create map...
  if (newvol->data) delete[] newvol->data;
  newvol->data = new float[newvol->xsize*newvol->ysize*newvol->zsize];
}

// Pad each side of the volmap's grid with zeros. Negative padding results
// in a trimming of the map
void pad(VolumetricData *vol, int padxm, int padxp, int padym, int padyp, int padzm, int padzp) {
  
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  double xdelta[3] = {(vol->xaxis[0] / (vol->xsize - 1)) , (vol->xaxis[1] / (vol->xsize - 1)) , (vol->xaxis[2] / (vol->xsize - 1))};
  double ydelta[3] = {(vol->yaxis[0] / (vol->ysize - 1)) , (vol->yaxis[1] / (vol->ysize - 1)) , (vol->yaxis[2] / (vol->ysize - 1))};
  double zdelta[3] = {(vol->zaxis[0] / (vol->zsize - 1)) , (vol->zaxis[1] / (vol->zsize - 1)) , (vol->zaxis[2] / (vol->zsize - 1))};
 

  int xsize_new = MAX(1, xsize + padxm + padxp);
  int ysize_new = MAX(1, ysize + padym + padyp);
  int zsize_new = MAX(1, zsize + padzm + padzp);
  
  int gridsize = xsize_new*ysize_new*zsize_new;
  float *data_new = new float[gridsize];
  memset(data_new, 0, gridsize*sizeof(float));

  int startx = MAX(0, padxm);
  int starty = MAX(0, padym);
  int startz = MAX(0, padzm);
  int endx = MIN(xsize_new, xsize+padxm);
  int endy = MIN(ysize_new, ysize+padym);
  int endz = MIN(zsize_new, zsize+padzm);

  int gx, gy, gz;
  for (gx=startx; gx<endx; gx++)
  for (gy=starty; gy<endy; gy++)
  for (gz=startz; gz<endz; gz++)
    data_new[gx + gy*xsize_new + gz*xsize_new*ysize_new] = vol->data[(gx-padxm) + (gy-padym)*xsize + (gz-padzm)*xsize*ysize];

  delete vol->data;
  vol->data = data_new;

  vol->xsize = xsize_new;
  vol->ysize = ysize_new;
  vol->zsize = zsize_new;
  

  vec_scaled_add(vol->xaxis, padxm+padxp, xdelta);
  vec_scaled_add(vol->yaxis, padym+padyp, ydelta);
  vec_scaled_add(vol->zaxis, padzm+padzp, zdelta);
  
  vec_scaled_add(vol->origin, -padxm, xdelta);
  vec_scaled_add(vol->origin, -padym, ydelta);
  vec_scaled_add(vol->origin, -padzm, zdelta);


}

void vol_com(VolumetricData *vol, float *com){
  double origin[3] = {0.0, 0.0, 0.0};
  double delta[3] = {0.0, 0.0, 0.0};
  origin[0] = vol->origin[0];
  origin[1] = vol->origin[1];
  origin[2] = vol->origin[2];
  delta[0] = (vol->xaxis[0] / (vol->xsize - 1)) + (vol->yaxis[0] / (vol->ysize - 1)) + (vol->zaxis[0] / (vol->zsize - 1));
  delta[1] = (vol->xaxis[1] / (vol->xsize - 1)) + (vol->yaxis[1] / (vol->ysize - 1)) + (vol->zaxis[1] / (vol->zsize - 1));
  delta[2] = (vol->xaxis[2] / (vol->xsize - 1)) + (vol->yaxis[2] / (vol->ysize - 1)) + (vol->zaxis[2] / (vol->zsize - 1));

  int xsize = vol->xsize;
  int ysize = vol->ysize;

  int gx,gy,gz;
  float ix,iy,iz;
  
  vec_zero(com);
  double mass = 0.0;

  for (int i = 0; i < vol->gridsize(); i++) {
    gz = i / (ysize*xsize);
    gy = (i % (ysize*xsize)) / xsize;
    gx = i % xsize;
    
    ix = origin[0] + (gx * delta[0]);
    iy = origin[1] + (gy * delta[1]);
    iz = origin[2] + (gz * delta[2]);
    float coord[3] = {ix,iy,iz};
    float m = vol->data[i];
    mass = mass+m;
    vec_scale(coord, m, coord);
    vec_add(com, com, coord); 
     
  }
  
  float scale = 1.0/mass;
  vec_scale(com, scale, com);
}

// Crop the map based on minmax values given in coordinate space. If
// the 'cropping box' exceeds the map boundaries, the map is padded
// with zeroes. 
void crop(VolumetricData *vol, double crop_minx, double crop_miny, double crop_minz, double crop_maxx, double crop_maxy, double crop_maxz) {

  double xdelta[3] = {(vol->xaxis[0] / (vol->xsize - 1)) , (vol->xaxis[1] / (vol->xsize - 1)) , (vol->xaxis[2] / (vol->xsize - 1))};
  double ydelta[3] = {(vol->yaxis[0] / (vol->ysize - 1)) , (vol->yaxis[1] / (vol->ysize - 1)) , (vol->yaxis[2] / (vol->ysize - 1))};
  double zdelta[3] = {(vol->zaxis[0] / (vol->zsize - 1)) , (vol->zaxis[1] / (vol->zsize - 1)) , (vol->zaxis[2] / (vol->zsize - 1))};

  // Right now, only works for orthogonal cells. look at cc_threaded in MDFF.C for example of calculating
  //correct delta for non-ortho
  int padxm = int((vol->origin[0] - crop_minx)/xdelta[0]);
  int padym = int((vol->origin[1] - crop_miny)/ydelta[1]);
  int padzm = int((vol->origin[2] - crop_minz)/zdelta[2]);

  int padxp = int((crop_maxx - vol->origin[0] - vol->xaxis[0])/xdelta[0]);
  int padyp = int((crop_maxy - vol->origin[1] - vol->yaxis[1])/ydelta[1]);
  int padzp = int((crop_maxz - vol->origin[2] - vol->zaxis[2])/zdelta[2]);

  pad(vol, padxm, padxp, padym, padyp, padzm, padzp);

}

void clamp(VolumetricData *vol, float min_value, float max_value) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  
  for (int i=0; i<xsize*ysize*zsize; i++) {
    if (vol->data[i] < min_value) vol->data[i] = min_value;
    else if (vol->data[i] > max_value) vol->data[i] = max_value;
  }
}



void scale_by(VolumetricData *vol, float ff) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 

  for (int i=0; i<xsize*ysize*zsize; i++)
    vol->data[i] *= ff;

}

void scalar_add(VolumetricData *vol, float ff) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 

  for (int i=0; i<xsize*ysize*zsize; i++)
    vol->data[i] += ff;

}

void fit_to_range(VolumetricData *vol, float min_value, float max_value) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 

  float min = vol->data[0];
  float max = vol->data[0];
  int i;
  for (i=1; i<xsize*ysize*zsize; i++) {
    if (vol->data[i] < min) min = vol->data[i];
    if (vol->data[i] > max) max = vol->data[i];
  }
  
  for (i=0; i<xsize*ysize*zsize; i++)
    vol->data[i] = min_value + (max_value - min_value)*(vol->data[i] - min)/(max - min);

}

// Downsample grid size by factor of 2
//Changed from volutil so it doesn't support PMF; unclear
//whether that is still important.
void downsample(VolumetricData *vol) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  int gx, gy, gz, j;
  
  int xsize_new = xsize/2;
  int ysize_new = ysize/2;
  int zsize_new = zsize/2;
  float *data_new = new float[xsize_new*ysize_new*zsize_new];
  
  int index_shift[8] = {0, 1, xsize, xsize+1, xsize*ysize, xsize*ysize + 1, xsize*ysize + xsize, xsize*ysize + xsize + 1};
  
  for (gx=0; gx<xsize_new; gx++)
  for (gy=0; gy<ysize_new; gy++)
  for (gz=0; gz<zsize_new; gz++) {
    int n_new = gx + gy*xsize_new + gz*xsize_new*ysize_new;
    int n = 2*(gx + gy*xsize + gz*xsize*ysize);
    double Z=0.;
    for (j=0; j<8; j++) Z += vol->data[n+index_shift[j]];
    data_new[n_new] = Z/8.;
  }
  
  xsize = xsize_new;
  ysize = ysize_new;
  zsize = zsize_new;
  //not stored in VolumetricData of VMD, only Volmap in volutil...
  //maybe we should add it because it is used so often.
  //vscale(xdelta, 2.);
  //vscale(ydelta, 2.);
  //vscale(zdelta, 2.);
  const float old_xaxis[3] = {(float)vol->xaxis[0], (float)vol->xaxis[1], (float)vol->xaxis[2]};
  const float old_yaxis[3] = {(float)vol->yaxis[0], (float)vol->yaxis[1], (float)vol->yaxis[2]};
  const float old_zaxis[3] = {(float)vol->zaxis[0], (float)vol->zaxis[1], (float)vol->zaxis[2]};
  
  float scaling_factor = 0.5*(xsize)/(xsize/2);
  vec_scale((float *)vol->xaxis, scaling_factor, old_xaxis);
  scaling_factor = 0.5*(ysize)/(ysize/2);
  vec_scale((float *)vol->yaxis, scaling_factor, old_yaxis);
  scaling_factor = 0.5*(zsize)/(zsize/2);
  vec_scale((float *)vol->zaxis, scaling_factor, old_zaxis);
  
//  vaddscaledto(origin, 0.25, xdelta);
//  vaddscaledto(origin, 0.25, ydelta);
//  vaddscaledto(origin, 0.25, zdelta);
      
  delete[] vol->data;
  vol->data = data_new;
}

// Supersample grid size by factor of 2
void supersample(VolumetricData *vol) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  
  int gx, gy, gz;
  int xsize_new = xsize*2-1;
  int ysize_new = ysize*2-1;
  int zsize_new = zsize*2-1;
  int xysize = xsize*ysize;
  int xysize_new = xsize_new*ysize_new;
  float *data_new = new float[xsize_new*ysize_new*zsize_new];
  
  // Copy map to the finer grid
  for (gx=0; gx<xsize; gx++)
    for (gy=0; gy<ysize; gy++)
      for (gz=0; gz<zsize; gz++)
        data_new[2*gx + 2*gy*xsize_new + 2*gz*xysize_new] = \
          vol->data[gx + gy*xsize + gz*xysize];

  // Perform cubic interpolation for the rest of the voxels

  // x direction
  for (gx=1; gx<xsize-2; gx++)
    for (gy=0; gy<ysize; gy++)
      for (gz=0; gz<zsize; gz++)
        data_new[2*gx+1 + 2*gy*xsize_new + 2*gz*xysize_new] = \
          cubic_interp(data_new[(2*gx-2) + 2*gy*xsize_new + 2*gz*xysize_new],
                       data_new[(2*gx)   + 2*gy*xsize_new + 2*gz*xysize_new],
                       data_new[(2*gx+2) + 2*gy*xsize_new + 2*gz*xysize_new],
                       data_new[(2*gx+4) + 2*gy*xsize_new + 2*gz*xysize_new],
                       0.5);
  // borders
  for (gy=0; gy<ysize; gy++)
    for (gz=0; gz<zsize; gz++) {
      // gx = 0
      data_new[1 + 2*gy*xsize_new + 2*gz*xysize_new] = \
        cubic_interp(data_new[0 + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[0 + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[2 + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[4 + 2*gy*xsize_new + 2*gz*xysize_new],
                     0.5);
      // gx = xsize-2
      data_new[2*(xsize-2)+1 + 2*gy*xsize_new + 2*gz*xysize_new] = \
        cubic_interp(data_new[2*(xsize-2)-2 + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[2*(xsize-2)   + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[2*(xsize-2)+2 + 2*gy*xsize_new + 2*gz*xysize_new],
                     data_new[2*(xsize-2)+2 + 2*gy*xsize_new + 2*gz*xysize_new],
                     0.5);
    }

  // y direction
  for (gx=0; gx<xsize_new; gx++)
    for (gy=1; gy<ysize-2; gy++)
      for (gz=0; gz<zsize; gz++)
        data_new[gx + (2*gy+1)*xsize_new + 2*gz*xysize_new] = \
          cubic_interp(data_new[gx + (2*gy-2)*xsize_new + 2*gz*xysize_new],
                       data_new[gx + (2*gy)*xsize_new   + 2*gz*xysize_new],
                       data_new[gx + (2*gy+2)*xsize_new + 2*gz*xysize_new],
                       data_new[gx + (2*gy+4)*xsize_new + 2*gz*xysize_new],
                       0.5);
  // borders
  for (gx=0; gx<xsize_new; gx++)
    for (gz=0; gz<zsize; gz++) {
      // gy = 0
      data_new[gx + 1*xsize_new + 2*gz*xysize_new] = \
        cubic_interp(data_new[gx + 0*xsize_new + 2*gz*xysize_new],
                     data_new[gx + 0*xsize_new + 2*gz*xysize_new],
                     data_new[gx + 2*xsize_new + 2*gz*xysize_new],
                     data_new[gx + 4*xsize_new + 2*gz*xysize_new],
                     0.5);
      // gy = ysize-2
      data_new[gx + (2*(ysize-2)+1)*xsize_new + 2*gz*xysize_new] = \
        cubic_interp(data_new[gx + (2*(ysize-2)-2)*xsize_new + 2*gz*xysize_new],
                     data_new[gx + 2*(ysize-2)*xsize_new     + 2*gz*xysize_new],
                     data_new[gx + (2*(ysize-2)+2)*xsize_new + 2*gz*xysize_new],
                     data_new[gx + (2*(ysize-2)+2)*xsize_new + 2*gz*xysize_new],
                     0.5);
    }

  // z direction
  for (gx=0; gx<xsize_new; gx++)
    for (gy=0; gy<ysize_new; gy++)
      for (gz=1; gz<zsize-2; gz++)
        data_new[gx + gy*xsize_new + (2*gz+1)*xysize_new] = \
          cubic_interp(data_new[gx + gy*xsize_new + (2*gz-2)*xysize_new],
                       data_new[gx + gy*xsize_new + (2*gz)*xysize_new],
                       data_new[gx + gy*xsize_new + (2*gz+2)*xysize_new],
                       data_new[gx + gy*xsize_new + (2*gz+4)*xysize_new],
                       0.5);
  // borders
  for (gx=0; gx<xsize_new; gx++)
    for (gy=0; gy<ysize_new; gy++) {
      // gz = 0
      data_new[gx + gy*xsize_new + 1*xysize_new] = \
        cubic_interp(data_new[gx + gy*xsize_new + 0*xysize_new],
                     data_new[gx + gy*xsize_new + 0*xysize_new],
                     data_new[gx + gy*xsize_new + 2*xysize_new],
                     data_new[gx + gy*xsize_new + 4*xysize_new],
                     0.5);
      // gz = zsize-2
      data_new[gx + gy*xsize_new + (2*(zsize-2)+1)*xysize_new] = \
        cubic_interp(data_new[gx + gy*xsize_new + (2*(zsize-2)-2)*xysize_new],
                     data_new[gx + gy*xsize_new + 2*(zsize-2)*xysize_new],
                     data_new[gx + gy*xsize_new + (2*(zsize-2)+2)*xysize_new],
                     data_new[gx + gy*xsize_new + (2*(zsize-2)+2)*xysize_new],
                     0.5);
  }


  vol->xsize = xsize_new;
  vol->ysize = ysize_new;
  vol->zsize = zsize_new;
  //vscale(xdelta, 0.5);
  //vscale(ydelta, 0.5);
  //vscale(zdelta, 0.5);

  delete[] vol->data;
  vol->data = data_new;
}

// Transform map to a sigma scale, so that isovalues in VMD correspond
// to number of sigmas above the mean
void sigma_scale(VolumetricData *vol) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  
  int size = xsize*ysize*zsize;
  double mean = 0.;
  int i;
  for (i=0; i<size; i++)
    mean += vol->data[i];
  mean /= size;

  double sigma = 0.;
  for (i=0; i<size; i++)
    sigma += (vol->data[i] - mean)*(vol->data[i] - mean);
  sigma /= size;
  sigma = sqrt(sigma);

  for (i=0; i<size; i++) {
    vol->data[i] -= mean;
    vol->data[i] /= sigma;
  }

}

// Makes a mask out of a map
// i.e. all values > 0 are set to 1
void binmask(VolumetricData *vol) {
  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  clamp(vol, 0., FLT_MAX);
  int i;
  for (i=0; i<xsize*ysize*zsize; i++) {
    if (vol->data[i] > 0) vol->data[i] = 1;
  }
}

//Gaussian blurring (as a 3D convolution), but the kernel can easily be changed to something else
void gauss3d(VolumetricData *vol, double sigma) {

  //Right now, only works if resolution is the same in all map dimensions
  double xdelta[3] = {(vol->xaxis[0] / (vol->xsize - 1)) , (vol->xaxis[1] / (vol->xsize - 1)) , (vol->xaxis[2] / (vol->xsize - 1))};

  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 

 
  // Pre-divide by sqrt(3) in x/y/z dimensions to get "sigma" in 3D
  sigma /= sqrt(3.);
  
  double delta = xdelta[0];
  int step = (int)(3.*sigma/delta); // size of gaussian convolution
  if (!step) return;
  
  // Build convolution kernel
  int convsize = 2*step+1;
  float *conv = new float[convsize*convsize*convsize];
  memset(conv, 0, convsize*convsize*convsize*sizeof(float));

/*
  // Pad the map if required
  if (flagsbits & USE_PADDING)
    pad(vol, convsize, convsize, convsize, convsize, convsize, convsize);
*/

  int gridsize = xsize*ysize*zsize;
  float *data_new = new float[gridsize];
  memset(data_new, 0, gridsize*sizeof(float));
  
  double r2, norm=0.;
  int cx, cy, cz; 
  for (cz=0; cz<convsize; cz++)
  for (cy=0; cy<convsize; cy++)
  for (cx=0; cx<convsize; cx++) {
    r2 = delta*delta*((cx-step)*(cx-step)+(cy-step)*(cy-step)+(cz-step)*(cz-step));
    conv[cx + cy*convsize + cz*convsize*convsize] = exp(-0.5*r2/(sigma*sigma)); 
    norm += conv[cx + cy*convsize + cz*convsize*convsize];
  }
  
  // Normalize...
  int n;
  for (n=0; n<convsize*convsize*convsize; n++) {
    conv[n] = conv[n]/norm;
  }
 
  // Apply convolution   
//  if (!ops->trivial())
  //  for (n=0; n<gridsize; n++) data[n] = ops->ConvertValue(data[n]);  
  
  int gx, gy, gz, hx, hy, hz; 
  for (gz=0; gz<zsize; gz++)
  for (gy=0; gy<ysize; gy++) 
  for (gx=0; gx<xsize; gx++)
  for (cz=0; cz<convsize; cz++)
  for (cy=0; cy<convsize; cy++)
  for (cx=0; cx<convsize; cx++) {
    hx=gx+cx-step;
    hy=gy+cy-step;
    hz=gz+cz-step;
    if (hx < 0 || hx >= xsize || hy < 0 || hy >= ysize || hz < 0 || hz >= zsize) {
      continue;
    }

    data_new[gx + gy*xsize + gz*xsize*ysize] += vol->data[hx + hy*xsize + hz*xsize*ysize]*conv[cx + cy*convsize + cz*convsize*convsize];  
  }
  
  //if (!ops->trivial())
  //  for (n=0; n<gridsize; n++) data_new[n] = ops->ConvertAverage(data_new[n]);  
  
  delete[] vol->data;
  vol->data = data_new;
}

// Fast Gaussian blur that takes advantage of the fact that the dimensions are separable.
void gauss1d(VolumetricData *vol, double sigma) {
  //Right now, only works if resolution is the same in all map dimensions
  double xdelta[3] = {(vol->xaxis[0] / (vol->xsize - 1)) , (vol->xaxis[1] / (vol->xsize - 1)) , (vol->xaxis[2] / (vol->xsize - 1))};

  int xsize = vol->xsize; 
  int ysize = vol->ysize; 
  int zsize = vol->zsize; 
  
  // Pre-divide by sqrt(3) in x/y/z dimensions to get "sigma" in 3D
  sigma /= sqrt(3.);
  
  double delta = xdelta[0];
  int step = (int)(3.*sigma/delta); // size of gaussian convolution
  if (!step) return;

  // Build convolution kernel
  int convsize = 2*step+1;
  float *conv = new float[convsize];
  memset(conv, 0, convsize*sizeof(float));
/*
  // Pad the map if required
  if (flagsbits & USE_PADDING)
    pad(convsize, convsize, convsize, convsize, convsize, convsize);
*/
  int gridsize = xsize*ysize*zsize;
  float *data_new = new float[gridsize];

  double r2, norm=0.;
  int c;
  for (c=0; c<convsize; c++) {
    r2 = delta*delta*(c-step)*(c-step);
    conv[c] = (float) exp(-0.5*r2/(sigma*sigma)); 
    norm += conv[c];
  }
  
  // Normalize...

  for (c=0; c<convsize; c++) {
    conv[c] = conv[c]/norm;
  }
/* 
  // Apply convolution   
  int n;
  if (!ops->trivial())
    for (n=0; n<gridsize; n++) data[n] = ops->ConvertValue(data[n]);  
*/
  memset(data_new, 0, gridsize*sizeof(float));
  
  int gx, gy, gz, hx, hy, hz; 
  for (gz=0; gz<zsize; gz++)
  for (gy=0; gy<ysize; gy++)
  for (gx=0; gx<xsize; gx++)
  for (c=0; c<convsize; c++) {
    hx=gx+c-step;
    hy=gy;
    hz=gz;
    if (hx < 0 || hx >= xsize) continue;
    data_new[gx + gy*xsize + gz*xsize*ysize] += vol->data[hx + hy*xsize + hz*xsize*ysize]*conv[c];
  }

  float *dataswap = vol->data;
  vol->data = data_new;
  data_new = dataswap;
  memset(data_new, 0, gridsize*sizeof(float));

  for (gz=0; gz<zsize; gz++)
  for (gy=0; gy<ysize; gy++)
  for (gx=0; gx<xsize; gx++)
  for (c=0; c<convsize; c++) {
    hx=gx;
    hy=gy+c-step;
    hz=gz;
    if (hy < 0 || hy >= ysize) continue;
    data_new[gx + gy*xsize + gz*xsize*ysize] += vol->data[hx + hy*xsize + hz*xsize*ysize]*conv[c];
  }
    
  dataswap = vol->data;
  vol->data = data_new;
  data_new = dataswap;
  memset(data_new, 0, gridsize*sizeof(float));
  
  for (gz=0; gz<zsize; gz++)
  for (gy=0; gy<ysize; gy++)
  for (gx=0; gx<xsize; gx++)
  for (c=0; c<convsize; c++) {
    hx=gx;
    hy=gy;
    hz=gz+c-step;
    if (hz < 0 || hz >= zsize) continue;
    data_new[gx + gy*xsize + gz*xsize*ysize] += vol->data[hx + hy*xsize + hz*xsize*ysize]*conv[c];
  }
/*
  if (!ops->trivial())
    for (n=0; n<gridsize; n++) data_new[n] = ops->ConvertAverage(data_new[n]);
  */
  delete[] vol->data;
  vol->data = data_new;
}

void add(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION) {

  int gx, gy, gz;

  // adding maps by spatial coords is slower than doing it directly, but
  // allows for precisely subtracting unaligned maps, and/or maps of
  // different resolutions

  if ( USE_UNION) {
    // UNION VERSION
    init_from_union(mapA, mapB, newvol);
  } else {
    // INTERSECTION VERSION
    init_from_intersection(mapA, mapB, newvol);
  }
  for (gx=0; gx<newvol->xsize; gx++)
    for (gy=0; gy<newvol->ysize; gy++)
      for (gz=0; gz<newvol->zsize; gz++) {
        float x, y, z;
        voxel_coord(gx, gy, gz, x, y, z, newvol);

        if (interp) newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_interpolate_from_coord(x,y,z) + \
          mapB->voxel_value_interpolate_from_coord(x,y,z);
        else newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_from_coord(x,y,z) + \
          mapB->voxel_value_from_coord(x,y,z);
      }

}

void subtract(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION) {

  int gx, gy, gz;

  // adding maps by spatial coords is slower than doing it directly, but
  // allows for precisely subtracting unaligned maps, and/or maps of
  // different resolutions

  if ( USE_UNION) {
    // UNION VERSION
    init_from_union(mapA, mapB, newvol);
  } else {
    // INTERSECTION VERSION
    init_from_intersection(mapA, mapB, newvol);
  }
  for (gx=0; gx<newvol->xsize; gx++)
    for (gy=0; gy<newvol->ysize; gy++)
      for (gz=0; gz<newvol->zsize; gz++) {
        float x, y, z;
        voxel_coord(gx, gy, gz, x, y, z, newvol);

        if (interp) newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_interpolate_from_coord(x,y,z) - \
          mapB->voxel_value_interpolate_from_coord(x,y,z);
        else newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_from_coord(x,y,z) - \
          mapB->voxel_value_from_coord(x,y,z);
      }

}

void multiply(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION) {

  int gx, gy, gz;

  // adding maps by spatial coords is slower than doing it directly, but
  // allows for precisely subtracting unaligned maps, and/or maps of
  // different resolutions

  if ( USE_UNION) {
    // UNION VERSION
    init_from_union(mapA, mapB, newvol);
  } else {
    // INTERSECTION VERSION
    init_from_intersection(mapA, mapB, newvol);
  }
  for (gx=0; gx<newvol->xsize; gx++)
    for (gy=0; gy<newvol->ysize; gy++)
      for (gz=0; gz<newvol->zsize; gz++) {
        float x, y, z;
        voxel_coord(gx, gy, gz, x, y, z, newvol);

        if (interp) newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_interpolate_from_coord(x,y,z) * \
          mapB->voxel_value_interpolate_from_coord(x,y,z);
        else newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          mapA->voxel_value_from_coord(x,y,z) * \
          mapB->voxel_value_from_coord(x,y,z);
      }

}

void average(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION) {

  int gx, gy, gz;

  // adding maps by spatial coords is slower than doing it directly, but
  // allows for precisely subtracting unaligned maps, and/or maps of
  // different resolutions

  if ( USE_UNION) {
    // UNION VERSION
    init_from_union(mapA, mapB, newvol);
  } else {
    // INTERSECTION VERSION
    init_from_intersection(mapA, mapB, newvol);
  }
  for (gx=0; gx<newvol->xsize; gx++)
    for (gy=0; gy<newvol->ysize; gy++)
      for (gz=0; gz<newvol->zsize; gz++) {
        float x, y, z;
        voxel_coord(gx, gy, gz, x, y, z, newvol);

        if (interp) newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          (mapA->voxel_value_interpolate_from_coord(x,y,z) + \
          mapB->voxel_value_interpolate_from_coord(x,y,z)) / 2;
        else newvol->data[gz*newvol->xsize*newvol->ysize + gy*newvol->xsize + gx] = \
          (mapA->voxel_value_from_coord(x,y,z) + \
          mapB->voxel_value_from_coord(x,y,z)) / 2;
      }

}

void vol_moveto(VolumetricData *vol, float *com, float *pos){
  float origin[3] = {0.0, 0.0, 0.0};
  origin[0] = (float)vol->origin[0];
  origin[1] = (float)vol->origin[1];
  origin[2] = (float)vol->origin[2];

  float transvector[3];
  vec_sub(transvector, pos, com);
  vec_add(origin, origin, transvector);   
  vol->origin[0] = origin[0];
  vol->origin[1] = origin[1];
  vol->origin[2] = origin[2];
}
/*
void vectrans(float *npoint, float *mat, double *vec){
  npoint[0]=vec[0]*mat[0]+vec[1]*mat[4]+vec[2]*mat[8];
  npoint[1]=vec[0]*mat[1]+vec[1]*mat[5]+vec[2]*mat[9];
  npoint[2]=vec[0]*mat[2]+vec[1]*mat[6]+vec[2]*mat[10];
}
*/
void vol_move(VolumetricData *vol,  float *mat){
  float origin[3] = {0.0, 0.0, 0.0};
  origin[0] = (float)vol->origin[0];
  origin[1] = (float)vol->origin[1];
  origin[2] = (float)vol->origin[2];
 
  float transvector[3] = {mat[12], mat[13], mat[14]};
 
  float npoint[3];
  
  //deal with origin transformation
  //vectrans  
  npoint[0]=origin[0]*mat[0]+origin[1]*mat[4]+origin[2]*mat[8];
  npoint[1]=origin[0]*mat[1]+origin[1]*mat[5]+origin[2]*mat[9];
  npoint[2]=origin[0]*mat[2]+origin[1]*mat[6]+origin[2]*mat[10];

  vec_add(origin, npoint, transvector);
  vol->origin[0] = origin[0]; 
  vol->origin[1] = origin[1]; 
  vol->origin[2] = origin[2]; 
     
  //deal with delta transformation
  double deltax[3] = {vol->xaxis[0],vol->xaxis[1],vol->xaxis[2]};
  double deltay[3] = {vol->yaxis[0],vol->yaxis[1],vol->yaxis[2]};
  double deltaz[3] = {vol->zaxis[0],vol->zaxis[1],vol->zaxis[2]};
 
  float npointx[3];
  float npointy[3];
  float npointz[3];
  vectrans(npointx, mat, deltax);
  vectrans(npointy, mat, deltay);
  vectrans(npointz, mat, deltaz);
  
  for (int i = 0; i<3; i++){
    vol->xaxis[i] = npointx[i];
    vol->yaxis[i] = npointy[i];
    vol->zaxis[i] = npointz[i];
  }

}

