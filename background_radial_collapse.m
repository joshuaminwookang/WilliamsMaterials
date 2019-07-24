% R-z collapses a given a z-leveled 3D pixel image.
% PRE: 1024 x 1024 x (___) image (most likely leveledImage) and the center of the image
% POST: returns the rz-collapse of the pixel intesities
%

function sampledIntensity = background_radial_collapse(image, scale, center)
    % declare variables and constants
    zDepth = size(image,3);
    pixCenter = center./scale(1,2) ;   % convert the center of circle to pixel coordinates
    
    % the maximum radius that we will be looking at
    rMax = floor(max([1024-pixCenter(1), 1024-pixCenter(2), pixCenter(1), pixCenter(2)]));
    
    % array to return; initialized generously in size just in case
    rz = zeros(floor(rMax*zDepth*0.7),2); 

    % array to return (final rz collapse)
    sampledIntensity = zeros(rMax,zDepth);
    temp = []; % temporary array to hold sampled intensities from 1~360 degrees

    % for each z, sample pixels on a line for 1~360 degrees to radially collapse image
    for z = 1:zDepth
        zSlice = permute(image(:,:,z),[2 1 3]); % the x-y plane to collpase
        for t = 1:360
            theta = t*pi/180;
            xEndPoints = rMax*cos(theta)*[0,1]+pixCenter(1);
            yEndPoints = rMax*sin(theta)*[0,1]+pixCenter(2);
            temp(t,:) = improfile(zSlice,xEndPoints,yEndPoints,rMax,'bicubic');
        end
        % take the average over all 360 degrees
        sampledIntensity(:,z) = mean(temp,'omitnan');
        temp = [];
    end
end
    
