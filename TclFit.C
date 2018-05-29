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



void moveby(AtomSel *sel, float *vect, MoleculeList *mlist, float *framepos){
  
  for (int i=sel->firstsel; i<=sel->lastsel; i++) {
    if (sel->on[i]) {
      vec_add(framepos + 3L*i, framepos + 3L*i, vect);
    }
  }
}

void move(AtomSel *sel, Matrix4 mat, MoleculeList *mlist, float *framepos){
  measure_move(sel, framepos, mat);
}

double calc_cc (AtomSel *sel, VolumetricData *volmapA, float resolution, MoleculeList *mlist, float *framepos) {
  
  float radscale;
  double gspacing = 0;
  double thresholddensity = 0.1;
  int verbose = 0;
  float return_cc = 0;

  radscale = .2*resolution;
  gspacing = 1.5*radscale;

  int quality = 0;
  if (resolution >= 9)
    quality = 0;
  else
    quality = 3;
  
  Molecule *mymol = mlist->mol_from_id(sel->molid());
  const float *radii = mymol->radius();
  
  int cuda_err = -1;
#if defined(VMDCUDA)
  VolumetricData **synthpp = NULL;
  VolumetricData **diffpp = NULL;
  VolumetricData **spatialccpp = NULL;

  if (getenv("VMDNOCUDA") == NULL) {

    cuda_err = vmd_cuda_compare_sel_refmap(sel, mlist, quality, 
                                  radscale, gspacing, 
                                  volmapA, synthpp, diffpp, spatialccpp, 
                                  &return_cc, thresholddensity, verbose);
  }
#endif

  // If CUDA failed, we use CPU fallback, and we have to prevent QuickSurf
  // from using the GPU either...
  if (cuda_err == -1) {
    const int force_cpu_calc=1;
    if (verbose)
      printf("Computing CC on CPUs...\n");

    if (gspacing == 0) {
      gspacing = 1.5*radscale;
    }

    QuickSurf *qs = new QuickSurf(force_cpu_calc);
    VolumetricData *volmap = NULL;
    volmap = qs->calc_density_map(sel, mymol, framepos, radii,
                                  quality, (float)radscale, (float)gspacing);
    double cc = 0.0;
    cc_threaded(volmap, volmapA, &cc, thresholddensity);
    return_cc += cc;
    delete qs;
  }
  return return_cc;
}

void do_rotate(int stride, float *com, AtomSel *sel, int amt, float *axis, MoleculeList *mlist, float *framepos){
  float move1[3];
  vec_scale(move1, -1.0, com);
  double amount = DEGTORAD(stride*amt);
  Matrix4 mat;
  mat.rotate_axis(axis, (float) amount);
  moveby(sel, move1, mlist, framepos);
  move(sel, mat, mlist, framepos);
  moveby(sel, com, mlist, framepos);

}

void rotate(int stride, int max_rot, float *com, float *returncc, float *bestpos, AtomSel *sel, MoleculeList *mlist, VolumetricData *volmapA, float resolution, float *origin, float *framepos) {
  
 // float *framepos = sel->coordinates(mlist);
 // float bestpos[sel->selected]; 
  //float best_cc = -1;
  for( int x = 0; x < max_rot; x++) {
    for( int y = 0; y < max_rot; y++) {
      for( int z = 0; z < max_rot; z++) {
        //x
        float axisx [3] = {1, 0, 0};
        do_rotate(stride, com, sel, x, axisx, mlist, framepos);
        //y
        float axisy [3] = {0, 1, 0};
        do_rotate(stride, com, sel, y, axisy, mlist, framepos);
        //z
        float axisz [3] = {0, 0, 1};
        do_rotate(stride, com, sel, z, axisz, mlist, framepos);
        
        float cc = calc_cc(sel, volmapA, resolution, mlist, framepos); 
        if (cc > *returncc) {
          *returncc = cc;
          
          for (int i=0; i<sel->selected*3L; i++) {
           // if (sel->on[i]) {
              bestpos[i] = framepos[i];
              framepos[i] = origin[i];
           // }
          }
        }
        for (int i=0; i<sel->selected*3L; i++) {
         // if (sel->on[i]) {
            framepos[i] = origin[i];
         // }
        }
       // framepos = origin;  
        
        //x
       // do_rotate(stride, com, sel, -x, axisx, mlist, framepos);
        //y
       // do_rotate(stride, com, sel, -y, axisy, mlist, framepos);
        //z
       // do_rotate(stride, com, sel, -z, axisz, mlist, framepos);

      }
    } 
  }
  //*returncc = best_cc;
  //return bestpos; 
}

void reset_origin(float *origin, float *newpos, AtomSel *sel) {
  for (int i=0; i<sel->num_atoms*3L; i++) {
   // if (sel->on[i]) {
      origin[i] = newpos[i];
   // }
  }
}

void calc_com(VolumetricData *vol, float *com){
  double origin[3] = {0.0, 0.0, 0.0};
  double delta[3] = {0.0, 0.0, 0.0};
  origin[0] = vol->origin[0];
  origin[1] = vol->origin[1];
  origin[2] = vol->origin[2];
  delta[0] = vol->xaxis[0] / (vol->xsize - 1);
  delta[1] = vol->yaxis[1] / (vol->ysize - 1);
  delta[2] = vol->zaxis[2] / (vol->zsize - 1);
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

void calc_moveto(VolumetricData *vol, float *com, float *pos){
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

void vectrans(float *npoint, float *mat, double *vec){
  npoint[0]=vec[0]*mat[0]+vec[1]*mat[4]+vec[2]*mat[8];
  npoint[1]=vec[0]*mat[1]+vec[1]*mat[5]+vec[2]*mat[9];
  npoint[2]=vec[0]*mat[2]+vec[1]*mat[6]+vec[2]*mat[10];
}

void calc_move(VolumetricData *vol,  float *mat){
  float origin[3] = {0.0, 0.0, 0.0};
  origin[0] = (float)vol->origin[0];
  origin[1] = (float)vol->origin[1];
  origin[2] = (float)vol->origin[2];
 
  float transvector[3] = {mat[12], mat[13], mat[14]};
  
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

int density_com(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
  int verbose = 0;
  if (argc < 3) {
     // "     options: --allframes (average over all frames)\n"
    Tcl_SetResult(interp, (char *) "usage: rigid "
      "com: [options]\n"
      "    options:  -i <input map> specifies new density filename to load.\n"
      "              -mol <molid> specifies an already loaded density's molid for use as target\n"
      "              -vol <volume id> specifies an already loaded density's volume id for use as target. Defaults to 0.\n",
//      "              --weight (weight density with atomic numbers)\n"
    //  "              -res <target resolution in Angstroms> (default 10.0)\n"
    //  "              -map <target resolution in Angstroms> (default 10.0)\n"
      TCL_STATIC);
    return TCL_ERROR;
  }


  int ret_val=0;
  int molid = -1;
  int volid = 0;
  const char *input_map = NULL;
  MoleculeList *mlist = app->moleculeList;
  
  for (int i=0; i < argc; i++) {
    char *opt = Tcl_GetStringFromObj(objv[i], NULL);
//    if (!strcmp(opt, "--weight")) {useweight = true;}
    if (!strcmp(opt, "-i")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No input map specified",NULL);
        return TCL_ERROR;
      }

      FileSpec spec;
      spec.waitfor = FileSpec::WAIT_BACK; // shouldn't this be waiting for all?
      input_map = Tcl_GetStringFromObj(objv[1+i], NULL);
      molid = app->molecule_new(input_map,0,1);
      //sel->molid()
      int ret_val = app->molecule_load(molid, input_map,app->guess_filetype(input_map),&spec);
      if (ret_val < 0) return TCL_ERROR;
    }


    if (!strcmp(opt, "-mol")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No molid specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &molid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n molid incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }

    if (!strcmp(opt, "-vol")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No volume id specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &volid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n volume id incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }


    if (!strcmp(opt, "-verbose") || (getenv("VMDMDFFVERBOSE") != NULL)) {
      verbose = 1;
    }
  }

  VolumetricData *volmapA = NULL;
  if (molid > -1) {
    Molecule *volmol = mlist->mol_from_id(molid);
    if (volmol == NULL) {
      Tcl_AppendResult(interp, "\n invalid molecule specified",NULL);
      return TCL_ERROR;
    }

    if (volmapA == NULL) 
      volmapA = volmol->modify_volume_data(volid);
  } else {
    Tcl_AppendResult(interp, "\n no target volume specified",NULL);
    return TCL_ERROR;
  }
  if (volmapA == NULL) {
    Tcl_AppendResult(interp, "\n no target volume correctly specified",NULL);
    return TCL_ERROR;
  }
  
  
  float com[3] = {0.0,0.0,0.0};
  calc_com(volmapA, com);  
 
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(com[0]));
  Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(com[1]));
  Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(com[2]));

  Tcl_SetObjResult(interp, tcl_result);
  return TCL_OK;
}

int density_move(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
  int verbose = 0;
  if (argc < 3) {
     // "     options: --allframes (average over all frames)\n"
    Tcl_SetResult(interp, (char *) "usage: rigid "
      "move: -mat <4x4 transform matrix to apply to density> [options]\n"
      "    options:  -i <input map> specifies new density filename to load.\n"
      "              -mol <molid> specifies an already loaded density's molid for use as target\n"
      "              -vol <volume id> specifies an already loaded density's volume id for use as target. Defaults to 0.\n"
      "              -o <filename> write moved density to file.\n",
//      "              --weight (weight density with atomic numbers)\n"
    //  "              -res <target resolution in Angstroms> (default 10.0)\n"
    //  "              -map <target resolution in Angstroms> (default 10.0)\n"
      TCL_STATIC);
    return TCL_ERROR;
  }


  int ret_val=0;
  int molid = -1;
  int volid = 0;
  float mat[16];
  const char *input_map = NULL;
  const char *outputmap = NULL;
  MoleculeList *mlist = app->moleculeList;
  
  for (int i=0; i < argc; i++) {
    char *opt = Tcl_GetStringFromObj(objv[i], NULL);
//    if (!strcmp(opt, "--weight")) {useweight = true;}
    if (!strcmp(opt, "-i")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No input map specified",NULL);
        return TCL_ERROR;
      }

      FileSpec spec;
      spec.waitfor = FileSpec::WAIT_BACK; // shouldn't this be waiting for all?
      input_map = Tcl_GetStringFromObj(objv[1+i], NULL);
      molid = app->molecule_new(input_map,0,1);
      //sel->molid()
      int ret_val = app->molecule_load(molid, input_map,app->guess_filetype(input_map),&spec);
      if (ret_val < 0) return TCL_ERROR;
    }


    if (!strcmp(opt, "-mol")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No molid specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &molid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n molid incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }

    if (!strcmp(opt, "-vol")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No volume id specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &volid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n volume id incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }
    
    if (!strcmp(opt, "-mat")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No matrix specified",NULL);
        return TCL_ERROR;
      } else if (tcl_get_matrix(Tcl_GetStringFromObj(objv[0],NULL), interp, objv[i+1], mat) != TCL_OK) {
          return TCL_ERROR;
      }
    }
    
    if (!strcmp(opt, "-o")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No output file specified",NULL);
        return TCL_ERROR;
      } else {
        outputmap = Tcl_GetStringFromObj(objv[i+1], NULL);
      }
    }

    if (!strcmp(opt, "-verbose") || (getenv("VMDMDFFVERBOSE") != NULL)) {
      verbose = 1;
    }
  }
  Molecule *volmol = NULL;
  VolumetricData *volmapA = NULL;
  if (molid > -1) {
    volmol = mlist->mol_from_id(molid);
    if (volmol == NULL) {
      Tcl_AppendResult(interp, "\n invalid molecule specified",NULL);
      return TCL_ERROR;
    }

    if (volmapA == NULL) 
      volmapA = volmol->modify_volume_data(volid);
  } else {
    Tcl_AppendResult(interp, "\n no target volume specified",NULL);
    return TCL_ERROR;
  }
  if (volmapA == NULL) {
    Tcl_AppendResult(interp, "\n no target volume correctly specified",NULL);
    return TCL_ERROR;
  }
  
 
  calc_move(volmapA, mat); 
  volmol->force_recalc(DrawMolItem::MOL_REGEN);

  if (outputmap != NULL) {
    volmap_write_dx_file(volmapA, outputmap);
  }
  return TCL_OK;

}

int density_moveto(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
  int verbose = 0;
  if (argc < 3) {
     // "     options: --allframes (average over all frames)\n"
    Tcl_SetResult(interp, (char *) "usage: rigid "
      "moveto: -pos <x,y,z coordinates to move com to> [options]\n"
      "    options:  -i <input map> specifies new density filename to load.\n"
      "              -mol <molid> specifies an already loaded density's molid for use as target\n"
      "              -vol <volume id> specifies an already loaded density's volume id for use as target. Defaults to 0.\n"
      "              -o <filename> write moved density to file.\n",
//      "              --weight (weight density with atomic numbers)\n"
    //  "              -res <target resolution in Angstroms> (default 10.0)\n"
    //  "              -map <target resolution in Angstroms> (default 10.0)\n"
      TCL_STATIC);
    return TCL_ERROR;
  }


  int ret_val=0;
  int molid = -1;
  int volid = 0;
  double pos[3] = {0.0, 0.0, 0.0};
  const char *input_map = NULL;
  const char *outputmap = NULL;
  MoleculeList *mlist = app->moleculeList;
  
  for (int i=0; i < argc; i++) {
    char *opt = Tcl_GetStringFromObj(objv[i], NULL);
//    if (!strcmp(opt, "--weight")) {useweight = true;}
    if (!strcmp(opt, "-i")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No input map specified",NULL);
        return TCL_ERROR;
      }

      FileSpec spec;
      spec.waitfor = FileSpec::WAIT_BACK; // shouldn't this be waiting for all?
      input_map = Tcl_GetStringFromObj(objv[1+i], NULL);
      molid = app->molecule_new(input_map,0,1);
      //sel->molid()
      int ret_val = app->molecule_load(molid, input_map,app->guess_filetype(input_map),&spec);
      if (ret_val < 0) return TCL_ERROR;
    }


    if (!strcmp(opt, "-mol")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No molid specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &molid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n molid incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }

    if (!strcmp(opt, "-vol")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No volume id specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &volid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n volume id incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }
    
    int num1;
    Tcl_Obj **vector;
    if (!strcmp(opt, "-pos")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No position coordinate specified",NULL);
        return TCL_ERROR;
      } else if (Tcl_ListObjGetElements(interp, objv[i+1], &num1, &vector) != TCL_OK) {
      return TCL_ERROR;
      }
    
      for (int i=0; i<num1; i++) {
        if (Tcl_GetDoubleFromObj(interp, vector[i], &pos[i]) != TCL_OK) {
          Tcl_SetResult(interp, (char *) "vecscale: non-numeric in vector", TCL_STATIC);
          return TCL_ERROR;
        }
      }
    
    }
    
    if (!strcmp(opt, "-o")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No output file specified",NULL);
        return TCL_ERROR;
      } else {
        outputmap = Tcl_GetStringFromObj(objv[i+1], NULL);
      }
    }

    if (!strcmp(opt, "-verbose") || (getenv("VMDMDFFVERBOSE") != NULL)) {
      verbose = 1;
    }
  }
  Molecule *volmol = NULL;
  VolumetricData *volmapA = NULL;
  if (molid > -1) {
    volmol = mlist->mol_from_id(molid);
    if (volmol == NULL) {
      Tcl_AppendResult(interp, "\n invalid molecule specified",NULL);
      return TCL_ERROR;
    }

    if (volmapA == NULL) 
      volmapA = volmol->modify_volume_data(volid);
  } else {
    Tcl_AppendResult(interp, "\n no target volume specified",NULL);
    return TCL_ERROR;
  }
  if (volmapA == NULL) {
    Tcl_AppendResult(interp, "\n no target volume correctly specified",NULL);
    return TCL_ERROR;
  }
  
 
  float com[3] = {0.0,0.0,0.0};
  float newpos[3] = {(float)pos[0], (float)pos[1], (float)pos[2]};
  calc_com(volmapA, com);  
  calc_moveto(volmapA, com, newpos);
  volmol->force_recalc(DrawMolItem::MOL_REGEN);

  if (outputmap != NULL) {
    volmap_write_dx_file(volmapA, outputmap);
  }
  return TCL_OK;

}

int fit(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp) {
  int verbose = 0;
  if (argc < 3) {
     // "     options: --allframes (average over all frames)\n"
    Tcl_SetResult(interp, (char *) "usage: rigid "
      "fit: <selection> -res <resolution of map in A> [options]\n"
      "    options:  -i <input map> specifies new target density filename to load.\n"
      "              -mol <molid> specifies an already loaded density's molid for use as target\n"
      "              -thresholddensity <x> (ignores voxels with values below x threshold)\n"
      "              -vol <volume id> specifies an already loaded density's volume id for use as target. Defaults to 0.\n",
//      "              --weight (weight density with atomic numbers)\n"
    //  "              -res <target resolution in Angstroms> (default 10.0)\n"
    //  "              -map <target resolution in Angstroms> (default 10.0)\n"
      TCL_STATIC);
    return TCL_ERROR;
  }

  //atom selection
  AtomSel *sel = NULL;
  sel = tcl_commands_get_sel(interp, Tcl_GetStringFromObj(objv[1],NULL));
  if (!sel) {
    Tcl_AppendResult(interp, "volmap: no atom selection.", NULL);
    return TCL_ERROR;
  }
  if (!sel->selected) {
    Tcl_AppendResult(interp, "volmap: no atoms selected.", NULL);
    return TCL_ERROR;
  }
  if (!app->molecule_valid_id(sel->molid())) {
    Tcl_AppendResult(interp, "invalide mol id.", NULL);
    return TCL_ERROR;
  }

  int ret_val=0;
  int molid = -1;
  int volid = 0;
  float radscale;
  double gspacing = 0;
  double resolution = 0;
  const char *input_map = NULL;
//  Tcl_GetDoubleFromObj(interp, objv[3], &resolution);
//  printf("Resolution %f\n", resolution);
//  const char *input_map = NULL;
  MoleculeList *mlist = app->moleculeList;
  Molecule *mymol = mlist->mol_from_id(sel->molid());
//  FileSpec spec;
//  spec.waitfor = FileSpec::WAIT_BACK; // shouldn't this be waiting for all?
//  input_map = Tcl_GetStringFromObj(objv[2], NULL);
//  molid = app->molecule_new(input_map,0,1);
  //sel->molid()
//  int ret_val = app->molecule_load(molid, input_map,app->guess_filetype(input_map),&spec);
//  if (ret_val < 0) return TCL_ERROR;
  
  for (int i=0; i < argc; i++) {
    char *opt = Tcl_GetStringFromObj(objv[i], NULL);
//    if (!strcmp(opt, "--weight")) {useweight = true;}
    if (!strcmp(opt, "-i")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No input map specified",NULL);
        return TCL_ERROR;
      }

      FileSpec spec;
      spec.waitfor = FileSpec::WAIT_BACK; // shouldn't this be waiting for all?
      input_map = Tcl_GetStringFromObj(objv[1+i], NULL);
      molid = app->molecule_new(input_map,0,1);
      //sel->molid()
      int ret_val = app->molecule_load(molid, input_map,app->guess_filetype(input_map),&spec);
      if (ret_val < 0) return TCL_ERROR;
    }

    if (!strcmp(opt, "-res")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No resolution specified",NULL);
        return TCL_ERROR;
      } else if (Tcl_GetDoubleFromObj(interp, objv[i+1], &resolution) != TCL_OK){ 
        Tcl_AppendResult(interp, "\n resolution incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }


    if (!strcmp(opt, "-mol")) {
      if (i == argc-1) {
        Tcl_AppendResult(interp, "No molid specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &molid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n molid incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }

    if (!strcmp(opt, "-vol")) {
      if (i == argc-1){
        Tcl_AppendResult(interp, "No volume id specified",NULL);
        return TCL_ERROR;
      } else if ( Tcl_GetIntFromObj(interp, objv[i+1], &volid) != TCL_OK) {
        Tcl_AppendResult(interp, "\n volume id incorrectly specified",NULL);
        return TCL_ERROR;
      }
    }


    if (!strcmp(opt, "-verbose") || (getenv("VMDMDFFVERBOSE") != NULL)) {
      verbose = 1;
    }
  }

  VolumetricData *volmapA = NULL;
  if (molid > -1) {
    Molecule *volmol = mlist->mol_from_id(molid);
    if (volmol == NULL) {
      Tcl_AppendResult(interp, "\n invalid molecule specified",NULL);
      return TCL_ERROR;
    }

    if (volmapA == NULL) 
      volmapA = volmol->modify_volume_data(volid);
  } else {
    Tcl_AppendResult(interp, "\n no target volume specified",NULL);
    return TCL_ERROR;
  }
  if (volmapA == NULL) {
    Tcl_AppendResult(interp, "\n no target volume correctly specified",NULL);
    return TCL_ERROR;
  }
/*
  // use quicksurf to compute simulated density map
  const float *framepos = sel->coordinates(app->moleculeList);
  const float *radii = mymol->radius();
  radscale = .2*resolution;

  if (gspacing == 0) {
    gspacing = 1.5*radscale;
  }

  int quality = 0;
  if (resolution >= 9)
    quality = 0;
  else
    quality = 3;

  if (verbose)
    printf("MDFF dens: radscale %f gspacing %f\n", radscale, gspacing);

  VolumetricData *synthvol=NULL;
  int cuda_err = -1;
#if defined(VMDCUDA)
  if (getenv("VMDNOCUDA") == NULL) {
    cuda_err = vmd_cuda_calc_density(sel, app->moleculeList, quality, radscale, gspacing, &synthvol, NULL, NULL, verbose);
    delete synthvol;
  }
#endif

  // If CUDA failed, we use CPU fallback, and we have to prevent QuickSurf
  // from using the GPU either...
  if (cuda_err == -1) {
    const int force_cpu_calc=1;
    QuickSurf *qs = new QuickSurf(force_cpu_calc);
    synthvol = qs->calc_density_map(sel, mymol, framepos, radii,
                                  quality, (float)radscale, (float)gspacing);
//    volmap_write_dx_file(volmap, outputmap);
    delete synthvol;
    delete qs;
  }
*/
  float *framepos = sel->coordinates(app->moleculeList);
  //compute center of mass 
  float com[3];
  // get atom masses
  const float *weight = mymol->mass();
  ret_val = measure_center(sel, framepos, weight, com);
  if (ret_val < 0) {
    Tcl_AppendResult(interp, "measure center failed",
         NULL);
    return TCL_ERROR;
  }
 
  float cc = -1;
  float *bestpos = new float [sel->num_atoms*3L];
  float *origin= new float [sel->num_atoms*3L];
  for (int k=0; k<sel->num_atoms*3L; k++) {
   // if (sel->on[i]) {
      origin[k] = framepos[k];
   // }
  }
  int stride = 24;
  int max_rot = 360/stride;
  rotate(stride, max_rot, com, &cc, bestpos, sel, mlist, volmapA, resolution, origin, framepos);
  reset_origin(origin, bestpos, sel); 

  int stride2 = 6;
  max_rot = stride/stride2;
  rotate(stride2, max_rot, com, &cc, bestpos, sel, mlist, volmapA, resolution, origin, framepos);
  rotate(-stride2, max_rot, com, &cc, bestpos, sel, mlist, volmapA, resolution, origin, framepos);
  reset_origin(origin, bestpos, sel); 
  
  int stride3 = 1;
  max_rot = stride2/stride3;
  rotate(stride3, max_rot, com, &cc, bestpos, sel, mlist, volmapA, resolution, origin, framepos);
  rotate(-stride3, max_rot, com, &cc, bestpos, sel, mlist, volmapA, resolution, origin, framepos);
  //reset_origin(origin, bestpos, sel); 
  
  //framepos = origin;
  for (int j=0; j<sel->num_atoms*3L; j++) {
   // if (sel->on[j]) {
      framepos[j] = bestpos[j];
   // }
  }
  //mymol->force_recalc(DrawMolItem::SEL_REGEN | DrawMolItem::COL_REGEN); 
  // notify molecule that coordinates changed.
  mymol->force_recalc(DrawMolItem::MOL_REGEN);
 
 // int frame = app->molecule_frame(sel->molid());
 // FileSpec speco;
 // speco.first = frame;                // write current frame only
 // speco.last = frame;                 // write current frame only
 // speco.stride = 1;                   // write all selected frames
 // speco.waitfor = FileSpec::WAIT_ALL; // wait for all frames to be written
 // speco.selection = sel->on;      // write only selected atoms
 // app->molecule_savetrajectory(sel->molid(), "fittedi.pdb", "pdb", &speco);
  printf("Best CC:%f\n", cc);
  return TCL_OK;
}

int obj_fit(ClientData cd, Tcl_Interp *interp, int argc,
                            Tcl_Obj * const objv[]){
  if (argc < 2) {
    Tcl_SetResult(interp,
    (char *) "Usage: rigid <command> [args...]\n"
      "Commands:\n"
      "fit      -- rigid body fitting\n"
      "com      -- get center of mass of density\n"
      "moveto   -- move density com to a specified coordinate\n"
      "move     -- apply specified 4x4 transformation matrix to density\n"
      ,
      TCL_STATIC);
    return TCL_ERROR;
  }
  char *argv1 = Tcl_GetStringFromObj(objv[1],NULL);

  VMDApp *app = (VMDApp *)cd;
  if (!strupncmp(argv1, "fit", CMDLEN))
    return fit(app, argc-1, objv+1, interp);
  if (!strupncmp(argv1, "com", CMDLEN))
    return density_com(app, argc-1, objv+1, interp);
  if (!strupncmp(argv1, "moveto", CMDLEN))
    return density_moveto(app, argc-1, objv+1, interp);
  if (!strupncmp(argv1, "move", CMDLEN))
    return density_move(app, argc-1, objv+1, interp);

  Tcl_SetResult(interp, (char *) "Type 'rigid' for summary of usage\n", TCL_VOLATILE);
  return TCL_OK;
}
