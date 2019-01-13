% This is a script to take stiffness data and find the slope of the contact
% data
%
% Justin Berman - January 4, 2017


% assumes you already have run the "measure_modulus" script, and have found
% the approximate beginning of the slope (over distance) and defined it as 
% a variable "start"

% end of where we are looking at the lope is 0.5 mm after the beginning
ending = start + 0.5;

% get all indecies of the data in dist where the distace in between the
% start and end
A = find(dist<ending & dist>start);

% get the slope going only 1 way
A = A(1:(length(A)/2));

% variables for the force and distance data in our range
mmN = force(A);
mmD = dist(A);

% Fit model to data.
[fitresult, gof] = fit(mmD, mmN, fittype( 'poly1' ));

% Plot fit with data.
plot(fitresult, mmD, mmN);
% Label axes
xlabel mmD
ylabel mmN
grid on

% get slope
slope = fitresult.p1;

% find E*
Estar = slope/(2*1.5e-3);
%find Young's Modulus
E = 3/4*Estar;

