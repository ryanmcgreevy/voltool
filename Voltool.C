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

