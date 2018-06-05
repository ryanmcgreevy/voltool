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
#include "Volutil.h"

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
//  printf("TRANSVEC: %f %f %f\n", mat[12],mat[13],mat[14]); 
  //const float dx = mat.mat[12];
  //const float dy = mat.mat[13];
  //const float dz = mat.mat[14];
 
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

