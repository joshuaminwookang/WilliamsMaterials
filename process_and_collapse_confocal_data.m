% Version - Summer 2018
%
% This is a simple script to process located confocal profile data, including:
%
% 1) pixel-to-micron conversion (if not already done)
% 2) zeroing / leveling the data by fitting subtracting the surrounding plane
% 3) r-z azimuthal collapse
% 4) make a nice plot of the collapsed data for this sphere
%
% requires the following .m files: 
% select_region.m, fit_poly22_to_points.m, radial_collapse_manual.m
%
% Kate Jensen - original assembled May 14, 2015 as a trial to give to Eva
%    revised and streamlined January 19, 2018 for Josh and Abdullah


%%%%%%%%%%%%
% NOTE: start with the data loaded in as 'r' (or change the variable names
% below)
%%%%%%%%%%%%

%% user-defined parameters 

% Need to change fileprefix!

%%%%%%%%%%%
% addpath /Users/jmk/Documents/Res_JensenLab/180703_Confocal_G81/180703_G81_Locating

% ask the user what to save these data as, e.g. G71sphere003
fileprefix = 'test';
%or if that's not working, you can just manually change things:
% fileprefix = 'sphere003';


%% (1) pixel-to-micrometer conversion
r = r_um;
%r = importdata('newTest.txt');
%(note that this overwrites r in the workspace)

%% (2) zeroing / leveling...


% zero/level the data by fitting a plane to the OUTER points, far from the
% indenting sphere
%
% In this case, we always put the sphere in the top-left-ish, so the
% bottom (high row number, large x) and right (high column number, large y)
% edges are where the data are flat. We'll use the points there to sample
% the flat plane.

% If that changes, you'll want to adjust the select_region inputs below.
% See >>help select_region for instructions on that function.

%only fit the edges of the image so that get a sampling of the
%undeformed plane around the particle
dist_from_edge = 10;
%for a symmetric box:     
%select_vector = ~select_region(r,'box',[dist_from_edge dist_from_edge -inf],[max(r(:,1))-dist_from_edge max(r(:,2))-dist_from_edge inf]);
%for only the right and bottom edges:
select_vector = select_region(r,'box',[0 0 -inf],[max(r(:,1))-dist_from_edge max(r(:,2))-dist_from_edge inf]);

%undeformed_plane_fit = fit_plane_to_points(r(~select_vector,:));
undeformed_plane_fit = fit_poly22_to_points(r(~select_vector,:)); %works better - takes out any curvature in surface

%generate zero surface to subtract at every point:
zero_surface_to_subtract = undeformed_plane_fit.p00 + undeformed_plane_fit.p10*r(:,1) + undeformed_plane_fit.p01*r(:,2)...
    + undeformed_plane_fit.p20*r(:,1).^2 + undeformed_plane_fit.p11*r(:,1).*r(:,2) + undeformed_plane_fit.p02*r(:,2).^2;

%subtract to level profile:
r_zeroed = [r(:,1:2) -zero_surface_to_subtract + r(:,3)];

%and save, put away
dlmwrite([fileprefix '_r_zeroed.txt'],r_zeroed,'\t')


%% (3) r-z azimuthal collapse

%now do the semi-manual collapse
% this requires manual inputs -- follow the instructions on the screen!
[center, rz] = radial_collapse_manual(r_zeroed);

%save the data for this sphere as text and .mat:
dlmwrite([fileprefix '_rz_collapsed.txt'],rz,'\t')
save([fileprefix '_rz_center_collapsed.mat'],'center','rz')


%% (4) make a nice plot of the collapsed data for this sphere

figure('Name',fileprefix)
plot(rz(:,1),rz(:,2),'.','DisplayName',['profile ' fileprefix])
hold all
%make the plot that looks nice:
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
axis equal
xlabel('r (\mu m) ')
ylabel('z (\mu m) ')
savefig([fileprefix '_collapsed_figure.fig'])

%% (5) collapse background fluorescence and plot
%[leveled,zero] = background_level_z(raw, scale, undeformed_plane_fit);
bg = background_rz_collapse(leveled, scale, center, zero);
%create plot
figure('Name',fileprefix)
plot(bg(:,1),bg(:,2),'.','DisplayName',['profile ' fileprefix])
hold all
plot(rz(:,1),rz(:,2),'.','DisplayName',['profile ' fileprefix])

%make the plot that looks nice:
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
%axis equal
axis ([0, 60, -15,5])
xlabel('r (\mu m) ')
ylabel('z (\mu m) ')
savefig([fileprefix '_bg_figure.fig'])


%% visualize-side view

for i=290:1024;
    imagesc(permute(leveled(i,:,:)>220,[3 2 1])); 
    hold all; 
    plot(rz(:,1)/0.14 + 330, rz(:,2)/0.225 + 61,'r.');
    hold off; 
    axis image; 
    axis xy; 
    title(i); 
    input(''); 
end