% want to integrate the confocal data I have, as solid of revolution, in
% order to compare the displaced volume under the sphere to the volume of
% the solid ridge around the sphere
%
% Kate Jensen - April 30, 2015
% tidied up to be a bit more general June 23, 2018

% this is a script, so assumes that the rz collapsed data is already loaded
% in with the following variables:
% R_d <-- a 2-column matrix listing the R and d
% radial_profiles <-- a cell array containing the collapsed profiles
%       This program runs through the list of profiles and pulls each out
%       as "rz" one at a time. If you only have one profile, you could also
%       just call it rz and skip to the integration below.


%load('180316_Confocal_G71S_Rd.mat');
%addpath /Users/jmk/Documents/Res_JensenLab/180316_Confocal_G71S/180316_G71S_Collapsed_Mat;
radial_profiles = cell(0);
%% load collapse_data %to use the example data
%load('G81_sphere003_rz_center_collapsed.mat');
radial_profiles{end+1} = rz;

%%
%run through all the profiles, and compute the volumes Delta_V_under and
%Delta_V_over

circle_color = [0.5 0.5 0.5]; %gray for fit circles

figure
for i = 1:length(radial_profiles)
    %% a bit of setup to display the mapped circle and the data we're working with:
    cla %clear current axes
    
    viscircles([0 R_d(i,1)-R_d(i,2)],R_d(i,1),'EdgeColor',circle_color,'LineWidth',4,'LineStyle','-');
    hold all
    rz = radial_profiles{i};
    plot(rz(:,1),rz(:,2),'.')
    
    
    axis equal
    axis([0 50 -20 10])
    set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold','PlotBoxAspectRatio',[4 3 1]); box on; grid on
    title(['R = ' num2str(R_d(i,1)) 'µm, d = ' num2str(R_d(i,2)) ' µm'])
    axis([0 60 -20 10]);
    xlabel('r  (µm)  ')
    ylabel('z  (µm)  ')
    
    %% prepare for integration:
    
    % in case there are outlier false-found points, remove everything above z=5µm
    rz(rz(:,2)>5,:)=[];
    
    %sort data by R so will integrate in the correct sequence; 
    %and make the very first entry be the sphere-fit-bottom at r=0 so the
    %integration doesn't randomly skip a chunk of the sphere volume 
    rz_to_integrate = [0 -R_d(i,2); sortrows(rz,1)];
    
    %now do a running average to figure out what's above and what's
    %below the undeformed plane
    clear z_smoothed
    l = 2; %how many points to smooth over
    for k = l+1:size(rz_to_integrate,1)-l
        z_smoothed(k,1) = mean(rz_to_integrate(k-l:k+l,2));
    end
    %find where the last negative value is
    I = find(z_smoothed>0,1,'first');
    
    %separate the two parts of the profile (under-over)
    rz_to_integrate_under = rz_to_integrate(1:I,:);
    rz_to_integrate_over = rz_to_integrate(I+1:end,:);
    
    %sanity check that I haven't messed up the profile
    hold all
  % plot(rz_to_integrate(1:end-l,1),z_smoothed,'--','LineWidth',2)
    plot(rz_to_integrate_under(:,1),rz_to_integrate_under(:,2),'o')
    plot(rz_to_integrate_over(:,1),rz_to_integrate_over(:,2),'s')
    
    %give me a chance to look at it and see that the over/under divisions make sense
    integrate_this(i,1) = input('integrate this profile? (1/0): ');
    
    %% if user says it's good to integrate, go for it!
    % (otherwise just move on to the next profile because there's something
    % wronge with this one -- or cancel out and fix the situation)
    
    if integrate_this(i,1)  
        %and run the integrations, as well as the analytic
        %calculations, and save the 
        Delta_V_under_analytic(i,1) = pi*R_d(i,2)^2/3*(3*R_d(i,1)-R_d(i,2)); %just from the geometry of the sphere (should be the same as the integrated data -- this is a check that the integration is working correctly)
        Delta_V_under(i,1) = integrate_data_solid_of_revolution(rz_to_integrate_under);
        Delta_V_over(i,1) = integrate_data_solid_of_revolution(rz_to_integrate_over);       
    end
    
end


%% finally, make some summary plots

% First, a plot to check the integration. The analytic and numerical
% calculations of the indented volume should be identical within
% experimental error, so this should be a straight line with slope 1:
figure('Name','Analytic vs. Numerical \DeltaV_{under}')
%figure(1); hold all
plot(Delta_V_under_analytic,-Delta_V_under,'.')
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on; % axis equal
xlabel('\DeltaV_{under} - analytic (µm^3)')
ylabel('\DeltaV_{under} - numerical (µm^3)')
axis([0 1e4 0 1e4])
title('should be perfect slope 1 line if integration ok')

% Next, the actual data -- what's the relationship between the excess
% (over) volume and the displaced volume?
%figure(2); hold all
figure('Name','\DeltaV_{over} vs. \DeltaV_{under}')
plot(-Delta_V_under,Delta_V_over,'.')
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on; % axis equal
xlabel('\DeltaV_{under} (µm^3)')
ylabel('\DeltaV_{over} (µm^3)')
axis([0 1e4 0 1e4])

% Finally, calculate what must be the liquid volume as the difference
% between the indented (under) volume and the ridge (over) volume, and plot
% that versus the indented volume, since that seems to be what mattered:
figure('Name','\DeltaV_{liquid} vs. \DeltaV_{under}')
loglog(-Delta_V_under,-Delta_V_under-Delta_V_over,'.')
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on; % axis equal
xlabel('\DeltaV_{under} (µm^3)')
ylabel('\DeltaV_{liquid} (µm^3)')
axis([1e2 1e5 1e1 1e5])


%% autosave what we've calculated in a temporary file
save 180316_G71S_displaced_volume.mat 

%%
figure('Name','\DeltaV_{liquid} vs. \DeltaV_{under}')
hold on;
loglog(-V_under,-V_under-V_over,'.');
loglog(-Delta_V_under,-Delta_V_under-Delta_V_over,'.');
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on;  % axis equal
set(gca,'xscale', 'log','yscale', 'log');
xlabel('\DeltaV_{under} (µm^3)')
ylabel('\DeltaV_{liquid} (µm^3)')
legend('Gelest 9:1', 'Gelest 7:1','show')








