% Function to refine the strain matrix and translation between two sets of
% tracked bead positions from the estimate given by stretch_est.
% When performing a time series of incremental stretches, it is also be
% possible to use the output of of stretch_refine at the previous time step
% as the new estimate.
% Written by Rob Style 14/11/2013
%
% [C,T]=stretch_refine(init_bead_pos,final_bead_pos,A,t,dist_between_pts,show_output)
%
% Inputs
%
% init_bead_pos
% These are the tracked bead positions in the initial state. Typically the
% output of the function centroids1
%
% final_bead_pos
% These are the tracked bead positions in the final state. Typically the
% output of the function centroids1
%
% A,t
% These are the outputs of stretch_est. A is an estimate of the strain
% matrix to get from initial to final bead positions. t is an estimate of
% the rigid body translation to get from initial to final bead positions
%
% dist_between_pts
% Optional parameter. The tracking algorithm compares the final bead
% position with an estimate of their final position using A,t. It tries to
% identify pairs of beads. If final bead positions are more than dist_between_pts from
% their estimated final position, the algorithm will not register that
% pair.
% If the estimate is good, or the points are dense, dist_between_pts should be small. If estimate is
% rough, then a larger dist_between_pts can be used. show_output (below)
% should be used to confirm accurate particle tracking.
% Default value is 2
%
% show_output
% Optional parameter. If this is set to 1, then a figure will display that
% allows you to check if the strain matrix is accurately calculated.
% Otherwise figure is suppressed.
% Default is 0.
%
% Outputs
% 
% C is the best fit 2x2 strain matrix that relates the initial points X=(x;y) to
% their final positions X'=(x';y'): X'=C.X+T
%
% T is the best fit translation vector that relates the initial points X=(x;y) to
% their final positions X'=(x';y'): X'=C.X+T

function [C,T]=stretch_refine(init_bead_pos,final_bead_pos,A,t,varargin)

if nargin==4
   show_output=0;
   dist_between_pts=2;
elseif nargin==5;
    dist_between_pts=varargin{1};
    show_output=0;
else
    dist_between_pts=varargin{1};
    show_output=varargin{2};
end

xyt=init_bead_pos;
xyt2=final_bead_pos;
xyt_est=xyt(:,1:2)*A'+repmat(t',length(xyt(:,1)),1);
trks=[xyt_est,ones(length(xyt_est(:,1)),1);xyt2(:,1:2),2*ones(length(xyt2(:,1)),1)];

param.good=2;
param.mem=0;
param.dim=2;
param.quiet=0;
tracks=track(trks,dist_between_pts,param);

xyt_trkd_est=tracks(1:2:end,1:2); % The estimated stretch points that are tracked to a final position point
xyt_trkd_final=tracks(2:2:end,1:2); % The final position points

% Now find the refining matrix that relates our guess x_g to the final
% points x_f: x_f=B x_g+t_2
BB=[xyt_trkd_est ones(length(xyt_trkd_est(:,1)),1)]\xyt_trkd_final;
B=BB(1:2,:)';
t2=BB(3,:)';

% Finally the best fit strain matrix and translation that goes from the
% initial points to the final points are C=B.A, T=B.t+t2
C=B*A;
T=B*t+t2;

if show_output==1
    figure
plot(xyt2(:,1),xyt2(:,2),'rx')
hold on
xyt_stretched=xyt(:,1:2)*C'+repmat(T',length(xyt(:,1)),1);
plot(xyt_stretched(:,1),xyt_stretched(:,2),'bo')
end

end