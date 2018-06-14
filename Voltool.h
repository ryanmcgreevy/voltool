#include "VolumetricData.h"
#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

void pad(VolumetricData *vol, int padxm, int padxp, int padym, int padyp, int padzm, int padzp);
void crop(VolumetricData *vol, double crop_minx, double crop_miny, double crop_minz, double crop_maxx, double crop_maxy, double crop_maxz);
void clamp(VolumetricData *vol, float min_value, float max_value);
void scale_by(VolumetricData *vol, float ff);
void scalar_add(VolumetricData *vol, float ff);
void fit_to_range(VolumetricData *vol, float min_value, float max_value);
void downsample(VolumetricData *vol);
void vol_com(VolumetricData *vol, float *com);
void vol_moveto(VolumetricData *vol, float *com, float *pos);
inline void vectrans(float *npoint, float *mat, double *vec){
  npoint[0]=vec[0]*mat[0]+vec[1]*mat[4]+vec[2]*mat[8];
  npoint[1]=vec[0]*mat[1]+vec[1]*mat[5]+vec[2]*mat[9];
  npoint[2]=vec[0]*mat[2]+vec[1]*mat[6]+vec[2]*mat[10];
}
void vol_move(VolumetricData *vol,  float *mat);
