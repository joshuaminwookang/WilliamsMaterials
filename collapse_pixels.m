
% Function to radially collapse background fluoroescence
% To be compared against rz-profile of the particles in order to find the 
% fluid volume on the surface

% PRE: center of image and the z-leveling vector has been found and given as a parameter
% POST: returns the rz-collapse of the pixel intesities
%
% Joshua Kang & Abdullah Nasir
% 2018. 11. 21
function bg = collapse_pixels(image, center, zero_surface_to_subtract, rz)
clear bg x y max_rad s

%find center in pixels 
x = round(center(1)/0.14);
y = round(center(2)/0.14);

%find max distance
max_rad = max([x+y, (1024-x)+y, (1024-y)+x, (1024-x)+(1024-y)]);

%intialize collapse, a 2D matrxi of r,z 
collapse(1:max_rad, 1:length(image(1,1,:))) = struct('avg',0,'count',0);

%for each z, iteratively compute average intensity based on input image
for k = 1: length(image(1,1,:)) 
    for i = 1 : length(image)
        for j = 1:length(image)
            rad = abs(i-x) + abs(j-y);
            if rad > 0
                sum = collapse(rad, k).avg * collapse(rad,k).count + image(i,j,k);
                collapse(rad, k).count = collapse(rad,k).count +1;
                collapse(rad, k).avg = sum/collapse(rad, k).count;
            end
        end
    end
end

%initialize bg, a 2D array of points (in um), where intensity is greater
% than the cutoff (@ 300)
bg = []
for i = 1:length(collapse(:,1)) 
    for j = 1:length(collapse(1,:))
        if collapse(i,j).avg > 300 
            bg = [bg; [i,j]];
        end
    end
end
bg = [0.14/1.3 0.225].* bg;
bg = [bg(:,1), bg(:,2) - mean(zero_surface_to_subtract)];

%save data
dlmwrite([fileprefix '_background_collapsed.txt'],bg,'\t')

