function pseudoGroupDifference(DTI_file_name, sphere_centers, L1_factor, ...
    ellipsoid_radii, noise_factor, FA_threshold, DTI_output_file_name)
%Alters the principle eigenvalues of NIFTI DTI data
%
% AUTHOR: Matthew Leming
% REQUIRED: vistasoft in your MATLAB path
%        (https://github.com/vistalab/vistasoft)
% DESCRIPTION: This function takes in the name of a NIFTI-1 Diffusion
%        Tensor Imaging file and increases the principle lambda values by a
%        certain factor. This function is meant to deliberately create
%        differences in data files in order to test difference detection
%        tools (such as TBSS).
% INPUTS:
%        - DTI_file_name - The path of the NIFTI-1 DTI file that this will
%          alter
%        - DTI_output_file_name - The name of the file to be written out
%        - radii_centers - A 3xN matlab matrix that represents the
%          coordinates of the spherical areas within the dataset that are
%          to be altered
%        - L1_factor - The number by which the principle eigenvalues
%          within each spherical area in the data ought to be multiplied
%          by
%        - noise_factor - the standard deviation of gaussian noise within
%          each sphere
%        - FA_threshold - if the FA of a certain voxel falls below this
%          value, it is unaffected. We don't want to do anything to CSF.

if exist('niftiWrite') ~= 2 || exist('niftiReadMatlab') ~= 2
    disp(' ')
    disp('   niftiWrite.m and/or niftiReadMatlab.m not found')
    disp('   Please get these from vistasoft and put those functions in')
    disp('   your MATLAB path. Run "addpath(genpath(/path/to/vistasoft))"')
    disp(' ')
    disp('   https://github.com/vistalab/vistasoft')
    disp(' ')
    return;
end    

if ~exist('DTI_file_name','var')
    disp('DTI_file_name is not set');
    return;
end
if ~exist('DTI_output_file_name','var')
     [pathstr,name,ext] = fileparts(DTI_file_name);
     DTI_output_file_name = sprintf('%s_altered.nii.gz',name);
end
if ~exist('sphere_centers','var')
     disp('sphere_centers not set')
     return;
elseif size(sphere_centers,1) ~= 3
    disp('sphere_centers is not of the correct dimensionality.')
    if size(sphere_centers,2) == 3 && size(size(sphere_centers),2) == 2
        disp('Transposing...')
        sphere_centers = sphere_centers';
    else
        return;
    end
end
if ~exist('L1_factor','var')
    L1_factor = 1.02;
    disp(sprintf('L1_factor not set. Setting to %f',L1_factor))
end
if ~exist('noise_factor','var')
    noise_factor = 0;
    disp(sprintf('noise_factor not set. Setting to %f',noise_factor))
end
if ~exist('ellipsoid_radii','var')
    ellipsoid_radii = 3;
    disp(sprintf('ellipsoid_radii not set. Setting to %d', ellipsoid_radii))
end
if ~exist('FA_threshold','var')
    FA_threshold = 0.01;
    disp(sprintf('FA_threshold not set. Setting to %f', FA_threshold))
end

nifti_data = niftiReadMatlab(DTI_file_name);
voxel_data = nifti_data.data;

% We multiple ellipsoid_data with every principle eigenvalue in the DTI
% dataset
ellipsoid_data = outputEllipsoidMultiNoise(...
    sphere_centers(1,:),sphere_centers(2,:),sphere_centers(3,:),...
    ellipsoid_radii,ellipsoid_radii,ellipsoid_radii,...
    nifti_data.dim,L1_factor, noise_factor);

% Preallocate these for speed
D11 = NaN;
D12 = NaN;
D13 = NaN;
D22 = NaN;
D23 = NaN;
D33 = NaN;
tensor_matrix = zeros([3 3]);
eigen_values = zeros(3,1);
FA = NaN;
mult_mat = NaN;


for x=1:nifti_data.dim(1)
    for y=1:nifti_data.dim(2)
        for z=1:nifti_data.dim(3)
            mult_factor = ellipsoid_data(x,y,z);
            if mult_factor ~= 1
                matrix_1_6 = voxel_data(x,y,z,1,:);
                D33 = matrix_1_6(1);
                D23 = matrix_1_6(2);
                D22 = matrix_1_6(3);
                D13 = matrix_1_6(4);
                D12 = matrix_1_6(5);
                D11 = matrix_1_6(6);
                tensor_matrix(:,:) =...
                   [D11 D12 D13 ;...
                    D12 D22 D23 ;...
                    D13 D23 D33];
                [vec,val] = eig(tensor_matrix);
                FA = sqrt(0.5) * sqrt((val(1,1) - val(2,2))^2 + ...
                                      (val(2,2) - val(3,3))^2 + ...
                                      (val(3,3) - val(1,1))^2)/ ...
                   sqrt(val(1,1)^2 + val(2,2)^2 + val(3,3)^2);
                if FA > FA_threshold
                    if val(1,1) > val(2,2)
                        if val(1,1) > val(3,3)
                            mult_mat = [ mult_factor 0 0 ; 0 1 0 ; 0 0 1 ];
                        else
                            mult_mat = [ 1 0 0 ; 0 1 0 ; 0 0 mult_factor ];
                        end
                    else
                        if val(2,2) > val(3,3)
                            mult_mat = [ 1 0 0 ; 0 mult_factor 0 ; 0 0 1 ];
                        else
                            mult_mat = [ 1 0 0 ; 0 1 0 ; 0 0 mult_factor ];
                        end
                    end
                    new_matrix_1_6 = vec * mult_mat * val * inv(vec);
                    voxel_data(x,y,z,1,:) = ...
                      [ new_matrix_1_6(3,3) ...
                        new_matrix_1_6(2,3) ...
                        new_matrix_1_6(2,2) ...
                        new_matrix_1_6(1,3) ...
                        new_matrix_1_6(1,2) ...
                        new_matrix_1_6(1,1) ];
                end
            end
        end
    end
end
nifti_data.data = voxel_data;
niftiWrite(nifti_data, DTI_output_file_name);