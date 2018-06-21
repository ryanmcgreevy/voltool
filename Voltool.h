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

inline void voxel_coord(int i, float &x, float &y, float &z, VolumetricData *vol){
  double xdelta[3] = {(vol->xaxis[0] / (vol->xsize - 1)), (vol->xaxis[1] / (vol->xsize - 1)), (vol->xaxis[2] / (vol->xsize - 1))};
  double ydelta[3] = {(vol->yaxis[0] / (vol->ysize - 1)), (vol->yaxis[1] / (vol->ysize - 1)), (vol->yaxis[2] / (vol->ysize - 1))};
  double zdelta[3] = {(vol->zaxis[0] / (vol->zsize - 1)), (vol->zaxis[1] / (vol->zsize - 1)), (vol->zaxis[2] / (vol->zsize - 1))};
  
  int xsize = vol->xsize;
  int ysize = vol->ysize;
   
  int gz = i / (ysize*xsize);
  int gy = (i % (ysize*xsize)) / xsize;
  int gx = i % xsize;
  
  x = vol->origin[0] + (gx * xdelta[0]) + (gy * xdelta[1]) + (gz * xdelta[2]);
  y = vol->origin[1] + (gx * ydelta[0]) + (gy * ydelta[1]) + (gz * ydelta[2]);
  z = vol->origin[2] + (gx * zdelta[0]) + (gy * zdelta[1]) + (gz * zdelta[2]);

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
