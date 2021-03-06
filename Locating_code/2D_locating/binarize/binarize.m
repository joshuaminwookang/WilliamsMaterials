function res = binarize(image_array,threshold,erosion,dilation)
% 
% NAME:
%               binarize
% PURPOSE:
%               binarizes image using built-in function imbinarize. Limited
%               capabilities to separate connected features into individual
%               entities using erosion and dilation of features.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = binarize( image_array, threshold, erosion, dilation )
% INPUTS:
%               image:  The two-dimensional array to be binarized (uint8).

%               threshold: brightness threshold for binarization

%               erosion: Number of Pixels the features will be eroded by.
%
%               dilation: Number of Pixels the features will be dilated by.
% OUTPUTS:
%               res:    binarized image.
% PROCEDURE:
%               Binarizes the image with a given threshold. The Particles
%               (white on black) can be eroded  and dilated to some degree
%               in order to separate particles that are close together. Any
%               holes in the brightness of the particles are filled to give
%               white single particles on black background.
% NOTES:

% MODIFICATION HISTORY: 
%                   Created by Hendrik Spanke Dec. 16


image_array = double(image_array);

BW = imbinarize(image_array,threshold); % Binary Image in BW

erode = strel('disk',erosion);  % create erosion mask for particles
BW = imerode(BW,erode);   % erode particles

dilate = strel('disk',dilation); % create dilation mask for particles
BW = imdilate(BW,dilate);    % dilate particles

BW = imfill(BW,'holes');  % fill in black pixels surrounded by white pixels.

res=BW;
end