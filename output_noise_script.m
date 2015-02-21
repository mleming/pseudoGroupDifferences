% A sample of how to run this, either minimalistically or not.

cy = [ 62 63 63 22 85 ];
cx = [ 65 52 78 71 56 ];
cz = [ 50 42 42 47 47 ];
ellipsoid_radii = 4;
FA_threshold = 0.01;
L1_factor = 1.02;
stddev = 0.005;
pseudoGroupDifference('DTI_QCMI_005_1_DTI_float.nii', [cx' cy' cz']',...
      L1_factor, stddev, ellipsoid_radii, FA_threshold,...
      'output_noise.nii.gz');
pseudoGroupDifference('DTI_QCMI_005_1_DTI_float.nii', [cx' cy' cz']')