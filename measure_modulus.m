% This is a script to load and analyze stiffness data from the Texture
% Analyzer.
%
% Kate Jensen and Justin Berman - December 7, 2017
% Improved by Jeremy - Aug 9, 2018


% assumes you already have the data loaded in as a Table, and then
% converted to an array via e.g.
% data = table2array(Gelest91);

a = importdata('SEAL_GEL8100_OvenCure24HR_1-3.tab');
data = a.data;
%%
type = 'SEAL_GEL8100_OvenCure24HR_1-3';
stiffness = [];
% the format of these data are 3 columns per test, consisting of 
% [Force (N), Distance (mm), Time (s)] (hopefully standard!)

%how many tests were there?
N = size(data,2); %<-- had better be an integer!
%%
for i = 1:3:N %(i directly indexes the Force column)
    force = data(:,i);
    dist = data(:,i+1);
    time = data(:,i+2);
    
    plot(dist,force,'.')
    %set starting point
    start = input('Where does it start?: ');
    fd_line_fit
    stiffness = cat(1, stiffness, [Estar,E]);
end

  E_sum = 0;
  Estar_sum =0;
for i = 1:N/3
    E_sum = E_sum + stiffness(i);
    Estar_sum = Estar_sum + stiffness(i,2);
end

    E_avg = E_sum/(N/3)
    Estar_avg = Estar_sum/(N/3)

save(['Bulk_stiffness_' type '.mat'], 'stiffness');

