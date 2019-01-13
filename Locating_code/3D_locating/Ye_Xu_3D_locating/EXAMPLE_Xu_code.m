% This is an example code for using the cntrd3 written by Ye Xu, cntrd3
% needs cntrd, pkfnd and bpass to run (these are in the 2d locating
% folder). The example tif stack that this code runs on is contained in
% Softliv_group_code/Locating_code/3D_tracking/Kilfoil_3D_locating

% The code works by finding the features in each frame in a tif stack,
% and finding the peak intensity of each of the features. It then tracks
% the same features where they appear in multiple frames, and looks at how
% the intensity changes from frame to frame. It then fits a gaussian or
% parabolic through this intensity to work out to sub-pixel level accuracy
% where the maximum intensity is for each feature.

% When the code analyses each frame, it first runs bpass on that frame.
% Then pkfnd, and finally cntrd. For tracking features between frames, it
% uses track.

% Note by Rob Style: The XY pixel accuracy can probably be improved by
% using Raghu radial symmetry tracking code (see the 2d locating folder)
% instead of pkfnd and cntrd.

clear all

% First we set all the parameters

% Find the tif file. If there multiple tifs, you can run a loop to to
% cntrd3 on each of the files that you find
im_folder='../Kilfoil_3D_locating/';
files=dir([im_folder,'*.tif']);
im_details=files(1);

% First we need to work out the size of the pixels in 3d. I would recommend
% doing this be opening the stack in imagej and seeing the diameter of the
% individual blobs, and how many frames in a row they appear:
% szx is a typical diameter of a blob
szx=5;

% szz is a typical number of frames in which a significant amount of
% feature appears
szz=11;

% Now we need to work out the threshold for the pkfnd function that is
% hidden within cntrd3. To do this, perform a bpass on a frame from the
% stack that contain some in focus features, and look at what a typical
% pixel value is inside the features.

frame_no=37;

img=imread([im_folder,im_details.name],frame_no);
im_filter=bpass(img,1,szx);

imagesc(im_filter)

%% Setting threshold and other parameters

% From the image, we see that a reasonable threshold for pkfnd (the image intensity of stuff
% outside the features) might be about 100. If blobs in the image are
% overlapping, then adjust szx to be smaller.
th=100;

% There are also a range of optional parameters:
%  param: a structure containing a few parameters
%       param.mem
%       param.quiet
%       param.mxdisp: maximum displacement of a particle from one slice to
%        the adjacent slice, used as maxdisp parameter in track function.
%        The default value is 1.
%       param.z1: the number of the first image in the stack you want to
%        analyze
%       param.z2: the number of the last image in the stack you want to
%        analyze
%       param.method: the method used to determine the 3d centroid.
%        "normal" uses the weighted average based on the intensity of the
%        particles in each slice.  "gaussian" fits the intensity curve with
%        a gaussion and find peaks in x, y, and z coordinates. I believe
%        gaussian is a better option
param.mxdisp = 1;
param.method = 'gaussian';
%% Run the cntrd3 and plot up the results
im_details.name=[im_folder,im_details.name];
pks = cntrd3(im_details,th,szx,szz,param);
pks(end,:)=[];
figure
plot3(pks(:,1),pks(:,2),pks(:,3),'.')

% NB please note that the Kilfoil 3d code is significantly faster than this
% code...