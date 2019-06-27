
% Given a 3D flourescent-dyed image, this function levels the background pixel 
% intensities in the z-direction
% PRE: 1024 x 1024 x (___) image (most likely, raw_residual) and a poly22 fit of the surface  
% POST: returns z-leveled 3D image
%

function [leveledImage, z_zero] = background_level_z(image, scale, plane_fit)

clear leveledImage X Y zShifts;

% for each (i,j) pixel coordinate, (1) convert to um (2) use plane_fit to find z_shift
% once z_shift is found, rescale by z' = z - z_shift
% POST: 3D matrix (image) leveled in the z-direction based on the plane_fit found from particles.


% pixel locations scaled to be actual distance (in um)
[X,Y] = meshgrid((1:size(image, 1))*scale(1), (1:size(image, 2))*scale(2));

% size of z_shifts for (x,y) coordinates 
zShifts = round((plane_fit.p00*ones(size(image,1)) + plane_fit.p10*X + plane_fit.p01*Y...
     + plane_fit.p11*(X.*Y) + plane_fit.p20*(X.*X) + plane_fit.p02*(Y.*Y))/scale(3));

z_zero = max(max(zShifts)); % find maximum z shift 
leveledImage = zeros(size(image,1), size(image,2), size(image,3)+2*z_zero); % intialize matrix to be returned

% fill up matrix to return with z-shifted values from original 3D image
for i = 1:size(image,1)
    for j = 1:size(image,2)
        zStack = image(i,j,:);
        leveledImage(i,j,z_zero - zShifts(i,j) + (1:size(image,3))) = zStack;
    end
end

