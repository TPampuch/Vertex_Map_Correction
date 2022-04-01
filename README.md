# Vertex_Map_Correction

## Why?
This R-package was created to create reliable maps based on the pairwise distance measurements, created with the Vertex Laser Geo 360 (Haglöf, Sweden) "compass" function. 
The Vertex provides the option to save GPS-coordinates, but in remote areas the distortion and thus the offset between coordinates can be very high. 
The solution is to measure "reference-trees" from more than one position (base) and use them to calculate and remove the offset. Manually, this can be 
a tiresome task, so here is an automated solution. 

## What?
The main function `VMC::VMC` needs the output of the Vertex Laser Geo 360 and the information on "base" (the position from which the tree was measured; 
which should be integers) and a column with the tree ID to find the overlapping trees.
For more information call `?VMC::VMC`
