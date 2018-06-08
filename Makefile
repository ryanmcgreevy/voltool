
# Makefile for vmd
# VMD for LINUXAMD64, version 1.9.4a19 (May 18, 2018)

.SILENT:

CD          = cd
COPY        = cp
COPYDIR     = cp -r
MAKEDIR     = mkdir -p
MOVE        = mv -f
DELETE      = rm -f
DELETEDIR   = rm -rf
ECHO        = echo
TAR         = tar
COMPRESS    = compress
GNUCOMPRESS = /usr/local/bin/gzip
LATEX       = latex

# makefile configuration
VPATH                   = .:../LINUXAMD64
SHELL                   = /bin/sh
.SUFFIXES: .C .c .cu ..c .i .o .ptx


# C++ source files
VMD_CCPP    =	OpenGLDisplayDevice.C \
		OpenGLPbufferDisplayDevice.C \
		OpenGLExtensions.C \
		OpenGLRenderer.C \
		OpenGLShader.C \
		OpenGLCache.C \
		colvaratoms.C \
		colvarbias.C \
		colvarbias_abf.C \
		colvarbias_alb.C \
		colvarbias_histogram.C \
		colvarbias_meta.C \
		colvarbias_restraint.C \
		colvar.C \
		colvarcomp.C \
		colvarcomp_angles.C \
		colvarcomp_coordnums.C \
		colvarcomp_distances.C \
		colvarcomp_protein.C \
		colvarcomp_rotations.C \
		colvardeps.C \
		colvargrid.C \
		colvarmodule.C \
		colvarparse.C \
		colvarproxy.C \
		colvarproxy_vmd.C \
		colvarscript.C \
		colvartypes.C \
		colvarvalue.C \
		IMDMgr.C \
		IMDSim.C \
		IMDSimThread.C \
		CmdIMD.C \
		imd.C \
		androidvmdstart.C \
		Animation.C \
		ArtDisplayDevice.C \
		AtomColor.C \
		AtomParser.C \
		AtomLexer.C \
		AtomRep.C \
		AtomSel.C \
		Axes.C \
		BaseMolecule.C \
		Benchmark.C \
		BondSearch.C \
		CmdAnimate.C \
		CmdColor.C \
		CmdDisplay.C \
		CmdLabel.C \
		CmdMaterial.C \
		CmdMenu.C \
		CmdMol.C \
		CmdRender.C \
		CmdTrans.C \
		CommandQueue.C \
		CoorPluginData.C \
		CUDAAccel.C \
		DisplayDevice.C \
		Displayable.C \
		DisplayRocker.C \
		DispCmds.C \
		DrawMolecule.C \
		DrawMolItem.C \
		DrawMolItem2.C \
		DrawMolItemRibbons.C \
		DrawMolItemMSMS.C \
		DrawMolItemNanoShaper.C \
		DrawMolItemRings.C \
		DrawMolItemOrbital.C \
		DrawMolItemQuickSurf.C \
		DrawMolItemSurface.C \
		DrawMolItemVolume.C \
		DrawForce.C \
		DrawRingsUtils.C \
		FastPBC.C \
		FileRenderList.C \
		FileRenderer.C \
		FPS.C \
		GaussianBlur.C \
		GelatoDisplayDevice.C \
		GeometryAngle.C \
		GeometryAtom.C \
		GeometryBond.C \
		GeometryDihedral.C \
		GeometryList.C \
		GeometryMol.C \
		GeometrySpring.C \
		Hershey.C \
		HMDMgr.C \
		Inform.C \
		ImageIO.C \
		Isosurface.C \
		JRegex.C \
		JString.C \
		macosxvmdstart.C \
		MaterialList.C \
		Matrix4.C \
		MayaDisplayDevice.C \
		MDFF.C \
		Measure.C \
		MeasureCluster.C \
		MeasurePBC.C \
		MeasureRDF.C \
		MeasureQCP.C \
		MeasureSurface.C \
		MeasureSymmetry.C \
		MeasureVolInterior.C \
		MobileInterface.C \
		MobileButtons.C \
		MobileTracker.C \
		Molecule.C \
		MoleculeList.C \
		MoleculeGraphics.C \
		MolFilePlugin.C \
		Mouse.C \
		MSMSInterface.C \
		NanoShaperInterface.C \
		Orbital.C \
		OrbitalJIT.C \
		PeriodicTable.C \
		P_JoystickTool.C \
		P_TugTool.C \
		P_GrabTool.C \
		P_PrintTool.C \
		P_PinchTool.C \
		P_UIVR.C \
		P_Buttons.C \
		P_Tracker.C \
		P_Tool.C \
		P_CmdTool.C \
		P_SensorConfig.C \
		ParseTree.C \
		PickList.C \
		PickModeAddBond.C \
		PickModeCenter.C \
		PickModeForce.C \
		PickModeList.C \
		PickModeMolLabel.C \
		PickModeMove.C \
		PickModeUser.C \
		PlainTextInterp.C \
		PluginMgr.C \
		POV3DisplayDevice.C \
		PSDisplayDevice.C \
		QMData.C \
		QMTimestep.C \
		QuickSurf.C \
		RadianceDisplayDevice.C \
		RayShadeDisplayDevice.C \
		R3dDisplayDevice.C \
		RenderManDisplayDevice.C \
		Scene.C \
		ScaleSpaceFilter.C \
		Segmentation.C \
		SnapshotDisplayDevice.C \
		Spaceball.C \
		SpaceballButtons.C \
		SpaceballTracker.C \
		SpatialSearch.C \
		SpringTool.C \
		Stage.C \
		STLDisplayDevice.C \
		Stride.C \
		Surf.C \
		SymbolTable.C \
		TachyonDisplayDevice.C \
		Timestep.C \
		UIObject.C \
		UIText.C \
		VMDApp.C \
		VMDCollab.C \
		VMDDir.C \
		VMDDisplayList.C \
		VMDMenu.C \
		VMDQuat.C \
		VMDTitle.C \
		VMDThreads.C \
		VolCPotential.C \
		VolMapCreate.C \
		VolMapCreateILS.C \
		VolumetricData.C \
		VolumeTexture.C \
		Volutil.C \
		VrmlDisplayDevice.C \
		Vrml2DisplayDevice.C \
		Watershed.C \
		WavefrontDisplayDevice.C \
		WKFThreads.C \
		WKFUtils.C \
		utilities.C \
		util_simd.C \
		vmd.C \
		vmdmain.C \
		X3DDisplayDevice.C \
		ColorInfo.C \
		TclCommands.C \
		TclFastPBC.C \
		TclMDFF.C \
		TclVoltool.C \
		TclMeasure.C \
		TclMolInfo.C \
		TclTextInterp.C \
		TclVec.C \
		TclGraphics.C \
		TclSegmentation.C \
		TclVolMap.C \
		cmd_animate.C \
		cmd_collab.C \
		cmd_color.C \
		cmd_display.C \
		cmd_imd.C \
		cmd_label.C \
		cmd_material.C \
		cmd_menu.C \
		cmd_mobile.C \
		cmd_mol.C \
		cmd_mouse.C \
		cmd_parallel.C \
		cmd_plugin.C \
		cmd_profile.C \
		cmd_render.C \
		cmd_spaceball.C \
		cmd_tool.C \
		cmd_trans.C \
		cmd_user.C \
		cmd_util.C \
		cmd_vmdbench.C \
		tcl_commands.C \
		VMDTkMenu.C

# C source files
VMD_CC      = 	hash.c \
		inthash.c \
		intstack.c \
		ptrstack.c \
		msmpot.c \
		msmpot_compute.c \
		msmpot_cubic.c \
		msmpot_setup.c \
		vmdsock.c \
		vmddlopen.c \
		pcre.c \
		fitrms.c

# CUDA source files
VMD_CU      = 	msmpot_cuda.cu \
		msmpot_cuda_latcut.cu \
		msmpot_cuda_shortrng.cu \
		CUDABench.cu \
		CUDAClearDevice.cu \
		CUDADispCmds.cu \
		CUDAFastPBC.cu \
		CUDAGaussianBlur.cu \
		CUDAMarchingCubes.cu \
		CUDAMDFF.cu \
		CUDAMeasureRDF.cu \
		CUDAMeasureQCP.cu \
		CUDAOrbital.cu \
		CUDAParPrefixOps.cu \
		CUDAQuickSurf.cu \
		CUDASegmentation.cu \
		CUDASpatialSearch.cu \
		CUDAUtil.cu \
		CUDAVolCPotential.cu \
		CUDAVolMapCreateILS.cu \
		CUDAWatershed.cu

# Header files
VMD_H       = 	OpenGLDisplayDevice.h \
		OpenGLPbufferDisplayDevice.h \
		OpenGLExtensions.h \
		OpenGLRenderer.h \
		OpenGLShader.h \
		OpenGLCache.h \
		colvar_UIestimator.h \
		colvaratoms.h \
		colvarbias.h \
		colvarbias_abf.h \
		colvarbias_alb.h \
		colvarbias_histogram.h \
		colvarbias_meta.h \
		colvarbias_restraint.h \
		colvarcomp.h \
		colvardeps.h \
		colvargrid.h \
		colvar.h \
		colvarmodule.h \
		colvarparse.h \
		colvarproxy.h \
		colvarproxy_vmd.h \
		colvarscript.h \
		colvartypes.h \
		colvarvalue.h \
		CUDAKernels.h \
		imd.h \
		IMDMgr.h \
		IMDSim.h \
		IMDSimThread.h \
		CmdIMD.h \
		Animation.h \
		ArtDisplayDevice.h \
		Atom.h \
		AtomColor.h \
		AtomParser.h \
		AtomRep.h \
		AtomSel.h \
		Axes.h \
		BaseMolecule.h \
		Benchmark.h \
		BondSearch.h \
		CmdAnimate.h \
		CmdColor.h \
		CmdDisplay.h \
		CmdLabel.h \
		CmdMaterial.h \
		CmdMenu.h \
		CmdMol.h \
		CmdRender.h \
		CmdTrans.h \
		Command.h \
		CommandQueue.h \
		CoorData.h \
		CUDAAccel.h \
		CoorPluginData.h \
		DepthSortObj.h \
		DispCmds.h \
		DisplayDevice.h \
		Displayable.h \
		DisplayRocker.h \
		DrawMolecule.h \
		DrawMolItem.h \
		DrawMolItemSolventPoints.data \
		DrawForce.h \
		GelatoDisplayDevice.h \
		FPS.h \
		FileRenderList.h \
		FileRenderer.h \
		Fragment.h \
		GeometryAngle.h \
		GeometryAtom.h \
		GeometryBond.h \
		GeometryDihedral.h \
		GeometryList.h \
		GeometryMol.h \
		GeometrySpring.h \
		Hershey.h \
		Inform.h \
		ImageIO.h \
		Isosurface.h \
		JRegex.h \
		JString.h \
		macosxvmdstart.h \
		MaterialList.h \
		Matrix4.h \
		MayaDisplayDevice.h \
		Measure.h \
		MeasureSymmetry.h \
		Molecule.h \
		MoleculeGraphics.h \
		MoleculeList.h \
		MolFilePlugin.h \
		Mouse.h \
		MSMSInterface.h \
		NameList.h \
		NanoShaperInterface.h \
		PeriodicTable.h \
		Orbital.h \
		P_JoystickTool.h \
		P_TugTool.h \
		P_PinchToo.h \
		P_GrabTool.h \
		P_PrintTool.h \
		P_Feedback.h \
		P_UIVR.h \
		P_Buttons.h \
		P_Tracker.h \
		P_CmdTool.h \
		P_SensorConfig.h \
		P_Tool.h \
		ParseTree.h \
		PickList.h \
		PickMode.h \
		PickModeAddBond.h \
		PickModeCenter.h \
		PickModeForce.h \
		PickModeList.h \
		PickModeMolLabel.h \
		PickModeMove.h \
		Pickable.h \
		PlainTextInterp.h \
		PluginMgr.h \
		PointerTool.h \
		POV3DisplayDevice.h \
		PSDisplayDevice.h \
		QMData.h \
		QMTimestep.h \
		RadianceDisplayDevice.h \
		RayShadeDisplayDevice.h \
		R3dDisplayDevice.h \
		ResizeArray.h \
		RenderManDisplayDevice.h \
		Residue.h \
		Scene.h \
		SnapshotDisplayDevice.h \
		SortableArray.h \
		Spaceball.h \
		SpaceballButtons.h \
		SpaceballTracker.h \
		SpatialSearch.h \
		SpringTool.h \
		Stack.h \
		Stage.h \
		STLDisplayDevice.h \
		Stride.h \
		Surf.h \
		SymbolTable.h \
		TachyonDisplayDevice.h \
		TextEvent.h \
		TextInterp.h \
		Timestep.h \
		UIObject.h \
		UIText.h \
		VMDApp.h \
		VMDDir.h \
		VMDDisplayList.h \
		VMDMenu.h \
		VMDQuat.h \
		VMDTitle.h \
		VMDThreads.h \
		VolCPotential.h \
		VolMapCreate.h \
		VolumetricData.h \
		VolumeTexture.h \
		Volutil.h \
		VrmlDisplayDevice.h \
		Vrml2DisplayDevice.h \
		Watershed.h \
		WavefrontDisplayDevice.h \
		X3DDisplayDevice.h \
		utilities.h \
		pcre.h \
		pcreinternal.h \
		pcretables.h \
		vmdsock.h \
		fitrms.h \
		TclCommands.h \
		TclTextInterp.h \
		tcl_commands.h \
		VMDTkMenu.h \
		plugin.h \
		molfile_plugin.h \
		libmolfile_plugin.h

# Header files
VMD_PTX     = 	

# yacc and lex files
VMD_YACC    = 	AtomParser.y

VMD_LEX     = 	AtomLexer.l

# Misc. data file
VMD_DATA    = 	.vmdsensors .vmdrc

VMD_EXTRA          = 	

VMD_OTHER_EXE      = 	../lib/stride/stride_LINUXAMD64 ../lib/surf/surf_LINUXAMD64 ../lib/tachyon/tachyon_LINUXAMD64

VMD_OTHER_NAMES    = 	stride_LINUXAMD64

VMD_MAIN_DIR_FILES = 	Announcement FEEDBACK LICENSE README configure

# Turn things into objects
VMD_OBJS    =   $(VMD_CCPP:.C=.o) $(VMD_CC:.c=.o) $(VMD_CU:.cu=.o)

INCDIRS     =         -I/usr/local/include/tcl -I../plugins/include -I../plugins/LINUXAMD64/molfile -I../lib/netcdf/include -I.

LIBS        = -lGL  -L/usr/X11R6/lib64 -lGL -lX11  -Wl,-rpath -Wl,$$ORIGIN/ -lcudart_static -lrt  -lXi -lpthread -ltcl8.5  -lmolfile_plugin -lnetcdf -lm -ldl $(VMDEXTRALIBS)

LIBDIRS     =     -L/usr/local/cuda-9.0/lib64    -L/usr/local/lib/tcl  -L../plugins/LINUXAMD64/molfile -L../lib/netcdf/lib_LINUXAMD64 

DEFINES     = -DVMDOPENGL -DVMDOPENGLPBUFFER -DVMDGLXPBUFFER   -DVMDCOLVARS -DVMDCUDA -DMSMPOT_CUDA -DVMDIMD -DVMDXINPUT -DVMDTHREADS -DWKFTHREADS -DUSEPOSIXTHREADS -D_REENTRANT -DVMDQUICKSURF -DVMDWITHCARBS -DVMDPOLYHEDRA -DVMDSURF -DVMDMSMS -DVMDNANOSHAPER -DVMDLATTICECUBES -DVMDTCL  -DVMDSTATICPLUGINS  

# compiler and compiler directives 
CC          = gcc
CFLAGS      = -m64 -Wall -Wno-unknown-pragmas -O6 -ffast-math -DARCH_LINUXAMD64 $(DEFINES) $(INCDIRS) 

CCPP	    = g++
CPPFLAGS    = -m64 -fno-for-scope -Wno-deprecated -Wall -Wno-unknown-pragmas -O6 -ffast-math  -DARCH_LINUXAMD64 $(DEFINES) $(INCDIRS)

NVCC        = nvcc
NVCCFLAGS   = -lineinfo --ptxas-options=-v -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_35 -gencode arch=compute_30,code=sm_37 -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=compute_70 -gencode arch=compute_70,code=sm_70 --ftz=true  --machine 64 -O3  -DARCH_LINUXAMD64 $(DEFINES) $(INCDIRS)

COMPILERC   = echo
RCFLAGS     = No resource compiler required on this platform.

DEPEND      = -MM
DEPENDFILE  = Makedata.depend

LOADLIBES   = $(LIBDIRS) $(LIBS) -rdynamic

LEX 	    = flex
YACC	    = yacc
YFLAGS      = -d

############################### 
# default rules 
###############################

.C.o: 
	$(ECHO) "Compiling " $< " --> " $*.o " ..."; \
	$(CCPP) $(CPPFLAGS) -c $< -o ../LINUXAMD64/$@

.c.o:
	$(ECHO) "Compiling " $< " --> " $*.o " ..."; \
	$(CC) $(CFLAGS) -c $< -o ../LINUXAMD64/$@

.cu.o:
	$(ECHO) "Compiling " $< " --> " $*.o " ..."; \
	$(NVCC) $(NVCCFLAGS) -c $< -o ../LINUXAMD64/$@

.cu.ptx:
	$(ECHO) "Compiling " $< " --> " $*.ptx " ..."; \
	$(NVCC) $(DEFINES) --use_fast_math -I/usr/local/encap/NVIDIA-OptiX-SDK-5.0.1-linux64/include -I/usr/local/cuda-9.0/include -gencode arch=compute_30,code=compute_30 -ptx $< -o ../LINUXAMD64/$@

.y.o:

.l.o:

########## Targets

### Source targets
all default:   vmd_LINUXAMD64

vmd_LINUXAMD64: y.tab.h $(VMD_OBJS) $(VMD_PTX)
	$(ECHO) "Linking " $@ "..."; \
	$(CD) ../LINUXAMD64 ; \
	$(CCPP) $(CPPFLAGS) -I../src -o $@ $(VMD_OBJS) $(LOADLIBES) ; 
	$(COMPILERC) $(RCFLAGS)

install:
	if [ ! -d "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ]; then \
		$(MAKEDIR) "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; \
	fi ; \
	if [ ! -d "/home/ryanmcgreevy/Dropbox/vmd/test" ]; then \
		$(MAKEDIR) "/home/ryanmcgreevy/Dropbox/vmd/test" ; \
	fi ; \
	if [ ! -d "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"/doc ]; then \
		$(MAKEDIR) "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"/doc; \
	fi
	-$(COPY) ../Announcement  "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	-$(COPY) ../README        "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	-$(COPY) ../LICENSE       "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	-$(COPY) ../doc/ug.pdf        "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"/doc
	if [ -f /home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd/vmd_LINUXAMD64 ]; then \
           $(MOVE) "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd/vmd_LINUXAMD64" "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd/OLD_vmd_LINUXAMD64" ; $(DELETE) "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd/OLD_vmd_LINUXAMD64" ; fi
	-$(COPY) ../LINUXAMD64/vmd_LINUXAMD64 "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	if [ -d "../lib/redistrib/lib_LINUXAMD64" ]; then \
		$(CD) ../lib/redistrib/lib_LINUXAMD64; $(TAR) -cf - ./* | \
		(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; $(TAR) -xf -) \
	fi ;
	-$(COPY) ../lib/stride/stride_LINUXAMD64 ../lib/surf/surf_LINUXAMD64 ../lib/tachyon/tachyon_LINUXAMD64 "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	-$(CD) ..; $(TAR) -cf - scripts | \
	(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; $(TAR) -xf -)
	-$(CD) ../lib; $(TAR) -cf - scripts | \
	(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; $(TAR) -xf -)
	-$(CD) ..; $(TAR) -cf - python | \
	(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"/scripts ; $(TAR) -xf -)
	-$(CD) ..; $(TAR) -cf - plugins | \
	(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; $(TAR) -xf -)
	-$(CD) ..; $(TAR) -cf - shaders | \
	(cd "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd" ; $(TAR) -xf -)
	if [ -f ../LINUXAMD64/OptiXShaders.ptx ]; then \
		$(COPY) ../LINUXAMD64/OptiXShaders.ptx "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd/shaders"; \
	fi; 
	-$(COPY) ../data/.vmdrc ../data/.vmdsensors ../data/vmd_completion.dat "/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"
	$(CD) ../bin ; \
	if [ -f run_vmd_tmp ]; then $(DELETE) run_vmd_tmp; fi ; \
	if [ ! -x "/bin/csh" ]; then \
		$(ECHO) "Info: /bin/csh shell not found, installing Bourne shell startup script instead" ; \
		$(ECHO) '#!/bin/sh' >> run_vmd_tmp ; \
		$(ECHO) 'defaultvmddir="/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"' >> run_vmd_tmp ; \
		$(ECHO) 'vmdbasename=vmd' >> run_vmd_tmp ; \
		cat vmd.sh >> run_vmd_tmp ; \
	else \
		$(ECHO) '#!/bin/csh' >> run_vmd_tmp ; \
		$(ECHO) 'set defaultvmddir="/home/ryanmcgreevy/Dropbox/vmd/test/lib/vmd"' >> run_vmd_tmp ; \
		$(ECHO) 'set vmdbasename=vmd' >> run_vmd_tmp ; \
		cat vmd.csh >> run_vmd_tmp ; \
	fi ; \
	chmod +x run_vmd_tmp ; \
	$(COPY) run_vmd_tmp "/home/ryanmcgreevy/Dropbox/vmd/test"/vmd ; \
	$(DELETE) run_vmd_tmp
	$(ECHO) Make sure "/home/ryanmcgreevy/Dropbox/vmd/test"/vmd is in your path.
	$(ECHO) "VMD installation complete.  Enjoy!"

##### remove most of the cruft
clean:
	$(CD) ../LINUXAMD64 ; \
		$(DELETE) *.ptx *.o *.lst *.warnings *.depend.old core

veryclean: clean
	$(CD) ../LINUXAMD64 ; \
	  $(DELETE) vmd_LINUXAMD64
	$(CD) ../src ; \
	  $(DELETE) *.ptx *.o *.lst *.a *~ core; \
	  $(DELETE) vmd_LINUXAMD64
	$(CD) ../doc ; \
	  $(DELETE) *~ core

# The '/usr/include' stuff is to reduce checking /usr/include dates
depend: y.tab.h
	if [ "$(DEPEND)" != "" ]; then \
	echo "Building Makefile dependencies"; \
	  $(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	  if [ -f $(DEPENDFILE) ]; then \
	    $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	  touch $(DEPENDFILE); \
	for i in ZZZ $(VMD_CCPP) ; do \
	  if [ "$$i" != "ZZZ" ]; then \
	        $(ECHO) checking dependencies for $$i ...; \
	        $(CCPP) $(DEPEND) $(CPPFLAGS) $$i |  \
			sed -e 's/\/usr\/include\/[^ ]*\..//g' \
			    -e 's/\.\.\/lib\/.*\/[^ ]*\..//g' | \
			grep -v '^ *\\$$' >> $(DEPENDFILE) ; \
	  fi; \
	done; \
	for i in ZZZ $(VMD_CC) ; do \
	  if [ "$$i" != "ZZZ" ]; then \
	        $(ECHO) checking dependencies for $$i ...; \
	        $(CC) $(DEPEND) $(CFLAGS) $$i |  \
			sed -e 's/\/usr\/include\/[^ ]*\..//g' \
			    -e 's/\.\.\/lib\/.*\/[^ ]*\..//g' | \
			grep -v '^ *\\$$' >> $(DEPENDFILE) ; \
	  fi; \
	done; \
	$(ECHO) ParseTree.o AtomLexer.o AtomParser.o: y.tab.h \
                >> $(DEPENDFILE); \
	$(COPY) $(DEPENDFILE) $(DEPENDFILE).LINUXAMD64 ; \
	else \
	  $(ECHO) "Cannot find dependencies; your compiler does not support dependency checking."; \
        fi



# to bootstrap without a Makedata.depend file, either
#   touch Makedata.depend
# or change the following line to 'sinclude'
include Makedata.depend

