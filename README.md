This was my personal git repository for the `voltool` plugin of [VMD](https://www.ks.uiuc.edu/Research/vmd/)
while I worked for the Theoretical and Computational Biophysics Group at the University of Illinois Urbana-Champaign.
`voltool` is a set of utility methods for manipulating the volumetric data (density maps) used by the `mdff` plugin of VMD.
The molecular dynamics flexible fitting (MDFF) method can be used to flexibly fit atomic structures into density maps (usually Cryo-EM, or
low-resolution X-ray crystallography for xMDFF). The method consists of adding external forces proportional to the gradient of the density map 
into a molecular dynamics (MD) simulation of the atomic structure.

This version of `voltool` represents a large overhaul of an older plugin to improve computational efficiency and add many new features including
a fast rigid-body fitting method.

I no longer work at UIUC or with current `voltool` or `mdff` development, so this repo is now archived. For current VMD related information please
see the [official site](https://www.ks.uiuc.edu/Research/vmd/)


# voltool
New rigid body fitting and volumetric data utilities for VMD.  
  
Usage: voltool \<command\> [args...]  
Commands:  
map operations using an atomic structure:  
fit          -- rigid body fitting  
cc           -- calculates the cross-correlation coefficient between a map and structure  
sim          -- creates a simulated map from an atomic structure  
operations on one map:  
com          -- get center of mass of density  
moveto       -- move density com to a specified coordinate  
move         -- apply specified 4x4 transformation matrix to density  
trim         -- trim edges of a density  
crop         -- crop density to values given in coordinate space  
clamp        -- clamp out of range voxel values  
smult        -- multiply every voxel by a scaling factor  
sadd         -- add a scaling factor to every voxel  
range        -- fit voxel values to a given range  
downsample   -- downsample by x2 (x8 total reduction)  
supersample  -- supersample by x2 (x8 total increase)  
sigma        -- transform map to sigma scale  
binmask      -- make a binary mask of the map  
smooth       -- 3D gaussian blur  
operations on two maps:  
add          -- add two maps together  
diff         -- subtract map2 from map1  
mult         -- multiply map1 and map2  
avg          -- average two input maps into one  
correlate    -- calculates the cross-correlation coefficient between two maps  
