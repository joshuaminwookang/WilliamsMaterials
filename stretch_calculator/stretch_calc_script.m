% Script to calculate the stretch between two beads images. Replace my
% example images with your own...


xyt=centroids1('Bottom_beads1_40xairnooptivar131024_155722L1Blue.tif',5,400);
xyt2=centroids1('Bottom_beads4_40xairnooptivar131024_164729L1Blue.tif',5,400);

[A,t]=stretch_est(xyt,xyt2); % Estimates strain and translation

[C,T]=stretch_refine(xyt,xyt2,A,t,2,1); % Refines strain and translation
% change the last input to 0 here   ^ if you don't want to see the figure
% at the end

% Finally calculate the principle stretches
E=C*C'; % Cauchy-Green tensor
lambda=sqrt(eig(E))-1; % Strains in the principle directions... i.e. principle stretches - 1