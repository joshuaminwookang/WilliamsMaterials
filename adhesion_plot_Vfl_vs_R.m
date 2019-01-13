% Version - Summer 2018
% by Joshua Kang 
% 
% Script to plot displaced fluid volume(V_fluid) 
% versus Sphere Radius x Young's modulus (R*E) 
% for soft adhesion confocal analysis 


% specify the Young's moduli of each silicone substrate 
% could be automated later, but for now ...

% load .mat file for each data 
% this part needs to be changed as we gather more data
load('180124_Confocal_Gelest91_Rd.mat')
G91 = R_d
load('180316_Confocal_G71S_Rd.mat')
G71S = R_d
load('180110_Confocal_G121L_Rd.mat')
G121L = R_d
load('180703_Confocal_Gelest81_Rd')
G81 = R_d
load('180316_Confocal_G91S_Rd.mat')
G91 = [G91;R_d]
load('180626_Confocal_DC_OC_Rd')
DC_OC = R_d
load('180627_Confocal_DC_RT_Rd')
DC_RT = R_d

load('displaced_volume_compiled.mat')
V_fl_G91 = cat(1,V_fl_G91, V_fl_G91S);
V_under_G91 = cat(1,V_under_G91, V_under_G91S);

%%
Vfl_vs_R(G91(:,1),V_fl_G91,G81(:,1),V_fl_G81,G121L(:,1),V_fl_G121,DC_OC(:,1),V_fl_DC_OC,DC_RT(:,1),V_fl_DC_RT)


%% 
%Vfl_vs_Vunder(V_under_G91,V_fl_G91,V_under_G121,V_fl_G121, V_under_DC_RT, V_fl_DC_RT, V_under_DC_OC, V_fl_DC_OC)
Vfl_vs_Vunder(V_under_G91,V_fl_G91,V_under_G121,V_fl_G121, V_under_G71, V_fl_G71,[],[])

%%
% Create figure
figure('Name','\DeltaV_{liquid} vs. \DeltaV_{under}');

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
scatter(V_under_G71,V_fl_G71,'DisplayName','Gelest 7:1','MarkerFaceColor','flat',...
    'MarkerEdgeColor','none');
scatter(V_under_G81,V_fl_G81,'DisplayName','Gelest 8:1','MarkerFaceColor','flat',...
    'MarkerEdgeColor','none');

% Create xlabel
xlabel('\DeltaV_{under} (µm^3)','FontWeight','bold');

% Create ylabel
ylabel('\DeltaV_{liquid} (µm^3)','FontWeight','bold');

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-5000 5000]);
box(axes1,'on');
grid(axes1,'on');
axis([0 6000 -1e5 1e5])
% Set the remaining axes properties
legend(axes1,'show');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',1);



