function [ ellipsoids ] = outputEllipsoidMultiNoise( cx, cy, cz,...
    rx, ry, rz,dimensions, avg_value, stddev)
% Writes a 3D array with ellipsoids throughout in designated spots
%
% AUTHOR: Matthew Leming
% DESCRIPTION: Used with pseudoGroupDifference.m in order to generate
%        arrays that can have factors with which to multiple each principle
%        lambda value in DTI data.
% INPUT:
%       - cx, cy, cz: one-dimensional arrays with the corresponding x-,y-,
%         and z-value centers for each ellipsoid in the array. Must be the
%         same dimensionality.
%       - rx, ry, rz: integers representing the radii of each ellipsoid.
%         Should probably be equal.
%       - dimensions: a 1x3 array with the dimensions of the matrix to be
%         created
%       - avg_value: the average value of the lambda multipliers within
%         each voxel; the rest of the values in the voxels are 1 (i.e.,
%         they have no effect on the diffusion.
%       - stddev: the standard deviation for the gaussian randomness
%         function that dictates how far from the mean each voxel within
%         the ellipsoids should be.

xrange = dimensions(1);
yrange = dimensions(2);
zrange = dimensions(3);

ellipsoids = ones([xrange yrange zrange]);

for i=1:size(cx,2)
    for x=1:xrange
        for y=1:yrange
            for z=1:zrange
                if ((x-cx(i))/rx)^2 + ...
                   ((y-cy(i))/ry)^2 + ...
                   ((z-cz(i))/rz)^2 < 1
                    ellipsoids(x,y,z) = normrnd(avg_value, stddev);
                end
            end
        end
    end
end

end

