% Code to show how to analyse image sequence with bwlabel and regionprops

% NAME:
%               bwlabel_regionprops_sequence_example
% PURPOSE:
%               Demonstrates how to use bwlabel and regionprops to find the
%               positions of particles in a sequence of images
% 
% CATEGORY:
%               Image Processing
%
% OUTPUTS:
%               pks:    A list of particle positions. First two columns are
%               particle positions. Third column is area of features as
%               determined with regionprops. Area allows for filtering of
%               features in a specific size-range.
%               frame number in the image sequence
% NOTES:
%               - First test out the parameters you need by running
%               bwlabel_regionprops_simplified_example.m on a typical image
%
% MODIFICATION HISTORY:
%               Written by Hendrik Spanke Dec 2016
%

%%

% Give the path to where the files are.
path = './Epi Scope Test/';

% Number of files. Let's assume that the files are called imageN.tif, where
% N are integers from 1:100
no_files=100;

% Give the three bpass parameters that you've done by running
% bwlabel_regionprops_simplified_example.m or bpass_gui.m on a typical image
a= [0 10 0];
% Similarly give the binarize parameters;
th=120;
erosion=2;
dilation=0;

% Prepare the variables that will hold the particle positions
p=[];
c=[];
pks=[];

% Now run a loop
for i=1:100
    im=imread([path,'image',num2str(i),'.tif']); % Read in the ith image
    im_pass=bpass(im,a(1),a(2),a(3)); % Do the bpass
    p=binarize(im_pass,th,sz,erosion,dilation); % Do the binarize
    pks=[pks; [c i*ones(length(c),1)]];
end

% HEY HENDRIK, WHEN YOU GET ROUND TO DOING THIS, COULD YOU USE WEIGHTED
% CENTROID INSTEAD OF CENTROID? I'VE UPDATED THE SIMPLE SCRIPT TO SHOW YOU
% WHAT I MEAN...