/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: VolumetricData.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.37 $	$Date: 2018/09/07 19:13:45 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Base class for storing volumetric data and associated gradient data
 *
 ***************************************************************************/

#ifndef VOLUMETRICDATA_H
#define VOLUMETRICDATA_H
#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))

/// Volumetric data class for potential maps, electron density maps, etc
class VolumetricData {
public:
  double origin[3];        ///< origin of volume (x=0, y=0, z=0 corner)
  double xaxis[3];         ///< direction and length for X axis (non-unit)
  double yaxis[3];         ///< direction and length for Y axis (non-unit)
  double zaxis[3];         ///< direction and length for Z axis (non-unit)
  int xsize, ysize, zsize; ///< number of samples along each axis
  char *name;              ///< human-readable volume dataset identifier
  float *data;             ///< raw data, total of xsize*ysize*zsize voxels
                           ///< stored x varying fastest, then y, then z.
  float *gradient;         ///< negated normalized volume gradient map
  float datamin, datamax;  ///< min and max data values 

  /// constructor
  VolumetricData(const char *name, const float *origin, 
                 const float *xaxis, const float *yaxis, const float *zaxis,
                 int xs, int ys, int zs, float *dataptr);

  VolumetricData(const char *name, const double *origin, 
                 const double *xaxis, const double *yaxis, const double *zaxis,
                 int xs, int ys, int zs, float *dataptr);

  /// destructor
  ~VolumetricData();

  /// return total number of gridpoints
  long gridsize() const { return long(xsize)*long(ysize)*long(zsize); }

  /// Sets data name to an internal copy of the provided string
  void set_name(const char* name);

  /// return cell side lengths
  void cell_lengths(float *xl, float *yl, float *zl) const;

  /// return cell axes
  void cell_axes(float *xax, float *yax, float *zax) const;

  /// return cell axes directions
  void cell_dirs(float *xax, float *yax, float *zax) const;

  /// return volumetric coordinate from cartesian coordinate
  void voxel_coord_from_cartesian_coord(const float *carcoord, float *voxcoord, int shiftflag) const;

  /// return index of the voxel nearest to a cartesian coordinate
  long voxel_index_from_coord(float xpos, float ypos, float zpos) const;

  /// return voxel at requested index, no safety checks
  inline float voxel_value(int x, int y, int z) const {
    return data[z*xsize*ysize + y*xsize + x];
  }

  /// return voxel, after safely clamping index to valid range
  float voxel_value_safe(int x, int y, int z) const;

  /// return interpolated value from 8 nearest neighbor voxels
  float voxel_value_interpolate(float xv, float yv, float zv) const;

  /// return voxel value based on cartesian coordinates
  float voxel_value_from_coord(float xpos, float ypos, float zpos) const;
  float voxel_value_interpolate_from_coord(float xpos, float ypos, float zpos) const;
  /// return voxel value based on cartesian coordinates. safe versions
  /// return zero if coordinates are outside the map.
  float voxel_value_from_coord_safe(float xpos, float ypos, float zpos) const;
  float voxel_value_interpolate_from_coord_safe(float xpos, float ypos, float zpos) const;


  /// (re)compute the volume gradient
  void compute_volume_gradient(void);

  /// provide the volume gradient
  void set_volume_gradient(float *gradient);

  /// return gradient at requested index, no safety checks
  void voxel_gradient_fast(int x, int y, int z, float *grad) const {
    long index = (z*xsize*ysize + y*xsize + x) * 3;
    grad[0] = gradient[index    ];
    grad[1] = gradient[index + 1];
    grad[2] = gradient[index + 2];
  }

  /// return gradient, after safely clamping index to valid range
  void voxel_gradient_safe(int x, int y, int z, float *grad) const;

  /// interpolate the gradient between the eight neighboring voxels
  void voxel_gradient_interpolate(const float *voxcoord, float *gradient) const;

  /// return voxel gradient based on cartesian coordinates
  void voxel_gradient_from_coord(const float *coord, float *gradient) const;
  void voxel_gradient_interpolate_from_coord(const float *coord, float *gradient) const;

  /// Cubic interpolation used by supersample
  inline float cubic_interp(float y0, float y1, float y2, float y3, float mu) {
    float mu2 = mu*mu;
    float a0 = y3 - y2 - y0 + y1;
    float a1 = y0 - y1 - a0;
    float a2 = y2 - y0;
    float a3 = y1;

    return (a0*mu*mu2+a1*mu2+a2*mu+a3);
  }
  
  /// get the cartesian coordinate of a voxel given its x,y,z indices
  inline void voxel_coord(int x, int y, int z, 
                          float &gx, float &gy, float &gz, 
                          VolumetricData *vol) {
    float xdelta[3], ydelta[3], zdelta[3];
    vol->cell_axes(xdelta, ydelta, zdelta);
    
    gx = vol->origin[0] + (x * xdelta[0]) + (y * ydelta[0]) + (z * zdelta[0]);
    gy = vol->origin[1] + (x * xdelta[1]) + (y * ydelta[1]) + (z * zdelta[1]);
    gz = vol->origin[2] + (x * xdelta[2]) + (y * ydelta[2]) + (z * zdelta[2]);
  }


  /// get the cartesian coordinate of a voxel given its 1D index 
  inline void voxel_coord(int i, float &x, float &y, float &z, 
                          VolumetricData *vol) {
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

  //
  //--unary ops--
  // 
  
  /// add or remove voxels in the given axis directions
  void pad(int padxm, int padxp, int padym, int padyp, int padzm, int padzp);
  
  /// crop a volumetric data to a given set of cartesian coordinates
  void crop(double crop_minx, double crop_miny, double crop_minz, double crop_maxx, double crop_maxy, double crop_maxz);

  /// clamp out of range voxel values
  void clamp(float min_value, float max_value);

  /// scales voxel data by given amount
  void scale_by(float ff);

  /// add scalar value to to all voxels
  void scalar_add(float ff);

  /// rescale voxel data to a given range
  void rescale_voxel_value_range(float min_value, float max_value);

  /// decimate/dowmnsample voxels by 2 in each dimension (x8 total reduction)
  void downsample();

  /// refine/supersample voxels by 2 in each dimension (x8 total increase)
  void supersample();

  /// Transform map to a sigma scale, so that isovalues in VMD correspond
  /// to number of sigmas above the mean
  void sigma_scale();

  /// Make a binary mask out of a map, i.e. all values > 0 are set to 1
  void binmask();
  
  ///Guassian Blur using algorithm from Segmentation
  void gaussian_blur(double sigma);
};


//
// Fast and loose accessor macros, don't use unless you have to
// 

/// fast but unsafe macro for querying volume gradients
#define VOXEL_GRADIENT_FAST_IDX(v, index, grad) \
  { (grad)[0] = v->gradient[index    ]; \
    (grad)[1] = v->gradient[index + 1]; \
    (grad)[2] = v->gradient[index + 2]; \
  }

/// fast but unsafe macro for querying volume gradients
#define VOXEL_GRADIENT_FAST(v, x, y, z, grad) \
  { long index = ((z)*v->xsize*v->ysize + (y)*v->xsize + (x)) * 3; \
    (grad)[0] = v->gradient[index    ]; \
    (grad)[1] = v->gradient[index + 1]; \
    (grad)[2] = v->gradient[index + 2]; \
  }

#endif // VOLUMETRICDATA_H
