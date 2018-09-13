/***************************************************************************
 *cr
 *cr            (C) Copyright 2007-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: TclMDFF.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.1 $      $Date: 2018/09/12 15:13:35 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Tcl bindings for MDFF functions
 *
 ***************************************************************************/

#ifndef TCLMDFF_H
#define TCLMDFF_H

#include "VMDApp.h"
#include <tcl.h>

int mdff_sim(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp);
int mdff_cc(VMDApp *app, int argc, Tcl_Obj * const objv[], Tcl_Interp *interp);

#endif
