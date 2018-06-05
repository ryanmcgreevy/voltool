void vol_com(VolumetricData *vol, float *com);
void vol_moveto(VolumetricData *vol, float *com, float *pos);
inline void vectrans(float *npoint, float *mat, double *vec){
  npoint[0]=vec[0]*mat[0]+vec[1]*mat[4]+vec[2]*mat[8];
  npoint[1]=vec[0]*mat[1]+vec[1]*mat[5]+vec[2]*mat[9];
  npoint[2]=vec[0]*mat[2]+vec[1]*mat[6]+vec[2]*mat[10];
}
void vol_move(VolumetricData *vol,  float *mat);
