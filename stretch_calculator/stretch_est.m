% Function to estimate the strain matrix and translation between two sets of
% tracked bead positions. A typical use would be to calculate the stretch
% and rotation of a sample of bottom beads of a TFM substrate when
% stretched.
% Written by Rob Style 14/11/2013
%
% stretch_est(init_bead_pos,final_bead_pos)
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
% Outputs
% 
% A is the best fit 2x2 strain matrix that relates the initial points X=(x;y) to
% their final positions X'=(x';y'): X'=A.X+t
%
% t is the best fit translation vector that relates the initial points X=(x;y) to
% their final positions X'=(x';y'): X'=A.X+t


function [A t]=stretch_est(init_bead_pos,final_bead_pos)

xyt=init_bead_pos;
xyt2=final_bead_pos;

figure
subplot(2,1,1)
plot(xyt(:,1),xyt(:,2),'rx')
title('Select 4 points, well spaced apart. Click a point in top screen')
subplot(2,1,2)
plot(xyt2(:,1),xyt2(:,2),'rx')

% Manually select 4 points and their final positions
for i=1:4
    subplot(2,1,1);
    [cee(i),ree(i)]=ginput(1);
    title(' ')
    subplot(2,1,2);
    title('Click same point in bottom figure')
    [cee2(i),ree2(i)]=ginput(1);
    title(' ')
    subplot(2,1,1)
    title('Click point in top figure')
end
title(' ')

% Find the points closest to where you clicked
for i=1:4
    [~,c]=min((xyt(:,1)-cee(i)).^2+(xyt(:,2)-ree(i)).^2);
    initpt(i,1:2)=xyt(c,1:2);
    [~,c]=min((xyt2(:,1)-cee2(i)).^2+(xyt2(:,2)-ree2(i)).^2);
    finalpt(i,1:2)=xyt2(c,1:2);
end

% Do least squares minimisation to get the estimated A and t
AA=[initpt ones(4,1)]\finalpt;
A=AA(1:2,:)';
t=AA(3,:)';

% Check that it looks reasonable.
figure
plot(xyt2(:,1),xyt2(:,2),'rx')
hold on
xyt_est=xyt(:,1:2)*A'+repmat(t',length(xyt(:,1)),1);
plot(xyt_est(:,1),xyt_est(:,2),'bo')
title('Points should line up on top of each other. If not, rerun')

end