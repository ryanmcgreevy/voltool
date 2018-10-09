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
