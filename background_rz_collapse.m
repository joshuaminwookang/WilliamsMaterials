
% Given a z-leveled 3D pixel image, r-z collapse 
% PRE: 1024 x 1024 x (___) image (most likely leveledImage) and the center of the image
% POST: returns the rz-collapse of the pixel intesities
%

function rz = background_rz_collapse(someImage, scale, center, z_zero)

pixCenter = center./scale(1,2) ;   % convert the center of circle to pixel coordinates

% the maximum radius that we will be looking at
rMax = floor(max([1024-pixCenter(1), 1024-pixCenter(2), pixCenter(1), pixCenter(2)]));
rz = zeros(rMax*size(someImage,3),3); % array to return 

% number of samples to take for each radius 
sampleNum = floor(rMax/sqrt(2))+1;

% resulting step size of the "sweeping angle" we will use to add up intensities 
% of pixels that lie on the circle of given radius
dTheta = (1/sampleNum) * pi/4;

% matrix representing degrees (think polar coordinates) of the pixel location to sample from 
D = repmat((0:pi/4:7/4*pi)',1,sampleNum) + repmat((0:dTheta:dTheta*(sampleNum-1)),8,1);
% matrix representing radii of the pixel locations to sample
R = repmat(reshape((1:rMax),[1,1,rMax]),size(D,1),size(D,2));

% row and columns of sampling positions in  pixel coordinates
row = floor(pixCenter(1) + R.*cos(D));
col = floor(pixCenter(2) + R.*sin(D));

% experimental code: try creating concentric circles with radii 1,2,3, ...
[radius, theta] = meshgrid(1:1:rMax, linspace(0,2*pi,100));
circlesX = pixCenter(1) + radius.*cos(theta);
circlesY = pixCenter(2) + radius.*sin(theta);


intensity = zeros(size(row,1)*size(row,2),1); %intialize pixel intensities array
counts = 1;
for z = 1:size(someImage,3)
    zSlice  = someImage(:,:,z);
    % find appropriate limit for circle radii
    for r = 1:rMax
        % fancy code to create a 2 x __ array of (i,j) coordinates of pixels to sample from;
        samplingPos = cat(3, row(:,:,r), col(:,:,r));
        samplingPos = permute(samplingPos, [2,3,1]);
        samplingPos = reshape(samplingPos, [size(samplingPos,1)*size(samplingPos,3), size(samplingPos,2)]);
        %for each row in our indexing array, record the pixel intensity 
        for l = 1:size(samplingPos,1)
            i = samplingPos(l,1); j = samplingPos(l,2);
            if (i<1 || j<1 || i>1024 || j>1024)
                intensity(l) = 0;
            else
                intensity(l) = zSlice(samplingPos(l,1),samplingPos(l,2));
            end
        end
        avg = sum(intensity )/length(intensity);
        if (avg >= 220)
            rz(counts,:) = [r*scale(1), (z-z_zero)*scale(3), avg]; 
            % record average intensity along with r,z values
            counts = counts+1;
        end
    end
end
