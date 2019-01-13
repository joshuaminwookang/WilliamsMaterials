% Version - Summer 2018
% by Joshua Kang 
% 
% Script to plot Sphere Indentation x Young's modulus(d*E) 
% versus Sphere Radius x Young's modulus (R*E) 
% for soft adhesion confocal analysis 

% Requires: dE_vs_RE.m

% specify the Young's moduli of each silicone substrate 
% could be automated later, but for now ...
stiffness = [5.11747;17.755653;11.08224252;11.49542619; 0.692025; 11.80515569; 0.2544]

% load .mat file for each data 
% this part needs to be changed as we gather more data
load('180124_Confocal_Gelest91_Rd.mat')
G91 = R_d
%load('180703_Confocal_Gelest81_Rd')
%G81 = R_d*stiffness(6,1)
load('180316_Confocal_G71S_Rd.mat')
G71S = R_d*stiffness(2,1)
load('180626_Confocal_DC_OC_Rd.mat')
DC_OC = R_d*stiffness(3,1)
load('180627_Confocal_DC_RT_Rd')
DC_RT = R_d*stiffness(4,1)
load('180110_Confocal_G121L_Rd')
G121L = R_d*stiffness(5,1)
load('180110_Confocal_GEL8100L_Rd')
GEL81 = R_d*stiffness(7,1)
load('180316_Confocal_G91S_Rd.mat')
G91 = [G91;R_d]*stiffness(1,1)




dE_vs_RE(G91(:,1),G91(:,2), G71S(:,1),G71S(:,2), G121L(:,1), G121L(:,2),[],[],[],[],[],[],[],[])
% plot dE vs RE (log-log)
%dE_vs_RE(G121L(:,1),G121L(:,2),G91(:,1),G91(:,2),G81(:,1),G81(:,2),G71S(:,1),G71S(:,2),DC_OC(:,1),DC_OC(:,2),DC_RT(:,1),DC_RT(:,2),GEL81(:,1),GEL81(:,2))
%% plot Fluid volume extraction vs RE