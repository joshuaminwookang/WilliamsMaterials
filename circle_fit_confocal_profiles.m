% 
% Version - Summer 2018
% A script to obtain sphere size and depth out of collapsed R-z profiles. 
% We'll then plot d vs R for the different data sets -- and
% see if we can fit them to get out a measurement of the surface stress
%
% - Prof Jensen 

% Edit by JMK

%% (by JMK) create final matrix of R vs d

% set appropriate name for data
fileprefix = 'Gelest81';
final = [];
addpath /Users/jmk/Documents/Res_JensenLab/180110_Confocal_G91L/180110_G91L_Colllapsed_RZ

% run this just once for your first rz profile data

%% first load in your collapsed rz profile as rz in the workspace

% Set current directory to where rz_collapsed_txt files are located
rz = dlmread('sphere008_rz_collapsed.txt') % change file name for each run!

%% determine (manually) the inner region (which will be fit to a circle)
figure('name','Confocal profiles')

plot(rz(:,1),rz(:,2),'.')
grid on; box on
set(gca,'LineWidth',1,'FontSize',16,'FontWeight','bold')

xlabel('Distance from indenter center (µm)','FontSize',18,'FontWeight','bold')
ylabel('Surface profile (µm)','FontSize',18,'FontWeight','bold')

xlim([0 40])

% this asks you to type in a number - what radial distance is definitely
% just the inner circle?
inner_region_limit = input('Inner region r-cutoff: ');


%% now fit a circle to the inner region of these to get the sphere position

%select only those points inside the desired R range, and mirror over the
%y-axis:
profile_to_fit = double(rz(~(rz(:,1)>inner_region_limit),:));
profile_to_fit = [[-profile_to_fit(:,1) profile_to_fit(:,2)]; profile_to_fit]; %now mirrored
% add these data dto the plot, highlighted:
hold all
plot(profile_to_fit(:,1), profile_to_fit(:,2),'.')


%try
    r0z0R_fit = circle_fit(profile_to_fit);
%catch
%    warning('error: data could not be fit')
%    r0z0R_fit = [0 0 0];
%end

R_d = [r0z0R_fit(:,3) (r0z0R_fit(:,3)-r0z0R_fit(:,2))];

display('Radius is: ')
R = r0z0R_fit(1,3)
display('indentation depth is: ')
d = (r0z0R_fit(:,3)-r0z0R_fit(:,2))

% I recommend you save these values into a two-column matrix R_d, e.g. using 
final= [final;R_d];
% where you make n the current sphere number. 
% Then you can easily plot d vs. R later on. 


% and add this result to the graph above so you can see if you are happy with it:
hold all
viscircles(r0z0R_fit(1,1:2),r0z0R_fit(1,3),'Color','k'); 
axis equal

%% Plot Final d vs R

% (1) Produce a d vs R plot based on data you just analyzed
% set appropriate name for figure
figure('name','180703 Gelest 8:1 R vs d')

scatter(final(:,1),final(:,2),'filled')
grid on; box on
set(gca,'LineWidth',1,'FontSize',12,'FontWeight','bold')

ylabel('Distance from Indenter Center, d (µm)','FontSize',14,'FontWeight','bold')
xlabel('Radius of Spehre, R (µm)','FontSize',14,'FontWeight','bold')

xlim([0 40])
ylim([0 20])

% (2) Save R vs d data for future use (compiled figure)
% change variable name !!!
R_d = final
save(['180703_Confocal_' fileprefix '_Rd.mat'], 'R_d');

