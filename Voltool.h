#include <stdio.h>
#include "VolumetricData.h"
#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

// Cubic interpolation used by supersample
inline float cubic_interp(float y0, float y1, float y2, float y3, float mu) {
  
  float mu2 = mu*mu;
  float a0 = y3 - y2 - y0 + y1;
  float a1 = y0 - y1 - a0;
  float a2 = y2 - y0;
  float a3 = y1;

  return (a0*mu*mu2+a1*mu2+a2*mu+a3);
}

inline void voxel_coord(int x, int y, int z, float &gx, float &gy, float &gz, VolumetricData *vol){
  float xdelta[3], ydelta[3], zdelta[3];
  vol->cell_axes(xdelta, ydelta, zdelta);
  
  gx = vol->origin[0] + (x * xdelta[0]) + (y * ydelta[0]) + (z * zdelta[0]);
  gy = vol->origin[1] + (x * xdelta[1]) + (y * ydelta[1]) + (z * zdelta[1]);
  gz = vol->origin[2] + (x * xdelta[2]) + (y * ydelta[2]) + (z * zdelta[2]);
 
}

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

//unary ops
void pad(VolumetricData *vol, int padxm, int padxp, int padym, int padyp, int padzm, int padzp);
void crop(VolumetricData *vol, double crop_minx, double crop_miny, double crop_minz, double crop_maxx, double crop_maxy, double crop_maxz);
void clamp(VolumetricData *vol, float min_value, float max_value);
void scale_by(VolumetricData *vol, float ff);
void scalar_add(VolumetricData *vol, float ff);
void fit_to_range(VolumetricData *vol, float min_value, float max_value);
void downsample(VolumetricData *vol);
void supersample(VolumetricData *vol);
void sigma_scale(VolumetricData *vol);
void binmask(VolumetricData *vol);
void gauss3d(VolumetricData *vol, double sigma);
void gauss1d(VolumetricData *vol, double sigma);


//binary ops
void add(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
void subtract(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
void multiply(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);
void average(VolumetricData *mapA, VolumetricData  *mapB, VolumetricData *newvol, bool interp, bool USE_UNION);

//Volumetric data utilities
void vol_com(VolumetricData *vol, float *com);
void vol_moveto(VolumetricData *vol, float *com, float *pos);
void vol_move(VolumetricData *vol,  float *mat);
void init_from_union(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol);
void init_from_intersection(VolumetricData *mapA, const VolumetricData *mapB, VolumetricData *newvol);
