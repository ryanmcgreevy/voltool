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

inline void vectrans(float *npoint, float *mat, double *vec){
  npoint[0]=vec[0]*mat[0]+vec[1]*mat[4]+vec[2]*mat[8];
  npoint[1]=vec[0]*mat[1]+vec[1]*mat[5]+vec[2]*mat[9];
  npoint[2]=vec[0]*mat[2]+vec[1]*mat[6]+vec[2]*mat[10];
}

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
void vol_com(VolumetricData *vol, float *com);
void vol_moveto(VolumetricData *vol, float *com, float *pos);
void vol_move(VolumetricData *vol,  float *mat);
