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
    % array to return; initialized generously in size just in caseß
    rz = zeros(floor(rMax*zDepth*0.7),2); 
    % a 360-by-rMax matrix where each column is 0,1,...,359 degrees for each radius
    sampledIntensity = zeros(rMax,zDepth);
    
    for z = 1:zDepth
        zSlice = image(:,:,z); % the x-y plane to collpase
        for theta = 1:pi/180:2*pi
            xEndPoints = rMax*cos(theta)*[0,1]+pixCenter(1)*[1,1];
            yEndPoints = rMax*sin(theta)*[0,1]+pixCenter(2)*[1,1];
            sampledIntensity(:,z) = sampledIntensity(:,z) + improfile(zSlice,xEndPoints,yEndPoints,rMax,'bicubic');
        end
        % take the average over all 360 degrees
        sampledIntensity(:,z) = sampledIntensity(:,z)/360;
    end

    % set cutoff intensity as 1.5*the intensity found at r=1, z=zDepth
    cutoff = sampledIntensity(1,zDepth);
    counter = 1; % how many 
    for z = 1:size(sampledIntensity,2)
        for r = 1:size(sampledIntensity,1)
            if(sampledIntensity(r,z) > cutoff) 
                rz(counter,:) = [r*scale(1),z*scale(3)];
            else
            end
        end
    end
end
