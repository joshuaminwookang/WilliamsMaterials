% This is a program to remove the effects of very uneven illumintation
% in the XY images.
%
% Use after the raw data has been read in (and inverted, if necessary), but
% before filtering the data.
%
% May not be appropriate for data with significant crystalline symmetry
%
% Written in particular for my collaboration with Peter Schall's group in
% Amsterdam, Summer 2013
%
% Kate Jensen - August 5, 2013


raw0 = raw; %also just for debugging

%when run, raw data already exists in the workspace as the 3D array 'raw'
%start by converting raw to be a double precision integer
raw = double(raw);
%measure the illumination variation by simply averaging up all the XY images
illumination_variation = mean(raw,3);
%smooth:
%smoothed_illumination_variation = smooth2(illumination_variation,5);
f = fspecial('disk');
smoothed_illumination_variation = imfilter(illumination_variation,f,'replicate');

%for debugging:
%figure; imagesc(smoothed_illumination_variation); axis image

%note the max value of the raw image; will restore this dynamic range:
raw_max_value = max(max(max(raw)));

%divide out the illumination variation - better to do as a for-loop, I
%think, than to build a new matrix the size of raw
for i = 1:size(raw,3)
    raw(:,:,i) = raw(:,:,i)./smoothed_illumination_variation;
end

%rescale raw to be 0-to-255
raw = raw.*255/max(max(max(raw)));

%for debugging:
figure; imagesc(smoothed_illumination_variation)
figure; imagesc(raw(:,:,1)); axis image
