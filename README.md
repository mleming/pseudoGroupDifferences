pseudoGroupDifferences
######################

This is a MATLAB script that alters the lambda values in a
number of spherical areas around a NIFTI diffusion tensor
imaging dataset. The purpose of this is to more accurately
test group population tools such as TBSS. To run this, one
simply needs the name of a NIFTI DTI file and a number of
coordinates of the spherical areas that one wants to alter.

    cy = [ 62 63 63 22 85 ];
    cx = [ 65 52 78 71 56 ];
    cz = [ 50 42 42 47 47 ];
    pseudoGroupDifference('data.nii', [cx' cy' cz']')

To see examples, check out output_noise_script.m

This requires that vistasoft be in your MATLAB path:

https://github.com/vistalab/vistasoft
