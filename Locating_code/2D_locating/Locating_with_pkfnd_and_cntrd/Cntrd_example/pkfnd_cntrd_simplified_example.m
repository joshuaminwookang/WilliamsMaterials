% Cleaned up code to demonstrate point tracking with pkfnd/cntrd

% NAME:
%               pkfnd_cntrd_simplified_example
% PURPOSE:
%               Demonstrates how to use pkfnd and cntrd to find the
%               positions of a bunch of particles
% 
% CATEGORY:
%               Image Processing
%
% OUTPUTS:
%               c:    A list of particle positions. First two columns are
%               particle positions. Other two columns are extra information
%               from cntrd (e.g. particle brightness)
% NOTES:
%               - Needs pkfnd gui to make more user friendly
%               - There should be an example image1.tif in this folder that
%               this example works with.
%               - For a sequence of images taken in the same lighting
%               conditions, you should do this on a typical image. Then use
%               the parameters on every image in the sequence in a loop...

% MODIFICATION HISTORY:
%               Written by Rob Style Dec 2016
%

%%

% Firstly work out the best parameters for the bpass
a=bpass_gui('before.tif');

% I found a =[0 10  0] to give quite nice looking results.
% NB try and keep the second parameter in bpass_gui as small as possible
% without shrinking the blobs down too much.
%%
im_pass=bpass(imread('image1.tif'),a(1),a(2),a(3));

% Ideally here, we would put in a gui to optimise the pkfnd parameters too.
% For the moment, pkfnd has two parameters that you need to input. The
% first is a threshold - a minimum intensity that might be a peak. We'll
% just estimate this by taking the average intensity of the picture

th=mean(mean(im_pass));

% The second parameter is a size that is slightly larger than then diameter
% of the average blobs in the image. This is about the same as 

sz=9;

p=pkfnd(im_pass,th,sz); % Gives a list with all the particle positions

% Display the results of pkfnd
imagesc(imread('image1.tif'))
hold on
plot(p(:,1),p(:,2),'rx')

%% Finally, we use cntrd to find the centres of the blobs to subpixel resolution.
% cntrd has one main input - a size which is a bit bigger than the average
% size of the blobs. It's recommended in the cntrd code that this is
% a(2)+2.
% So let's just try sz+2
% Is it worth having a gui here to optimise this? I think it's probably ok
% to use the results from a pkfnd gui and then either take sz or sz+2?

c=cntrd(im_pass,p,sz+2);

% Display the results

imagesc(imread('image1.tif'))
hold on
plot(c(:,1),c(:,2),'rx')
