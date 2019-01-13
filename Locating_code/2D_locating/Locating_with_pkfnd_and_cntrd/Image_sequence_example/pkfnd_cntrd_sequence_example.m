% Cleaned up code to show how to analyse image sequence with pkfnd/cntrd

% NAME:
%               pkfnd_cntrd_sequence_example
% PURPOSE:
%               Demonstrates how to use pkfnd and cntrd to find the
%               positions of particles in a sequence of images
% 
% CATEGORY:
%               Image Processing
%
% OUTPUTS:
%               pks:    A list of particle positions. First two columns are
%               particle positions. Third/fourth columns are extra information
%               from cntrd (e.g. particle brightness). Final column is the
%               frame number in the image sequence
% NOTES:
%               - First test out the parameters you need by running
%               pkfnd_cntrd_simplified_example.m on a typical image
%
% MODIFICATION HISTORY:
%               Written by Rob Style Dec 2016
%

%%

% Give the path to where the files are.
path = './Epi Scope Test/';

% Number of files. Let's assume that the files are called imageN.tif, where
% N are integers from 1:100
no_files=100;

% Give the three bpass parameters that you've done by running
% pkfnd_cntrd_simplified_example.m or bpass_gui.m on a typical image
a= [0 10 0];
% Similarly give the pkfnd parameters;
th=400;
sz=9;
% Ditto with cntrd parameters;
sz_cntrd=9;

% Prepare the variables that will hold the particle positions
p=[];
c=[];
pks=[];

% Now run a loop
for i=1:100
    im=imread([path,'image',num2str(i),'.tif']); % Read in the ith image
    im_pass=bpass(im,a(1),a(2),a(3)); % Do the bpass
    p=pkfnd(im_pass,th,sz);
    c=cntrd(im_pass,p,sz_cntrd);
    pks=[pks; [c i*ones(length(c),1)]];
end
