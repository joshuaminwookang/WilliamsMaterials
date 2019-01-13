% integrate data - solid of revolution
%
% code to take r-z profile data and integrate it numerically as a solid of
% revolution
%
% Kate Jensen - April 29, 2015

% when this becomes a function, rz will be an input. *the* input.

%assume there's a variable in the workspace called rz, which is the profile
%data. doesn't have to be sorted; just the data is fine. should be in true
%cylindrical coordinates, though, so r>=0 always. Actually, r(1) = 0, 
% z(1) = z(r=0) will usually be standard. Or r(1)= r0 such that z(r1) = 0

%sort the data so that consecutive r values are adjacent

function V = integrate_data_solid_of_revolution(rz)

rz_sorted = sortrows(rz,1);

V = pi/2 * sum( (rz_sorted(1:end-1,2)+rz_sorted(2:end,2)).*(rz_sorted(2:end,1).^2-rz_sorted(1:end-1,1).^2) );

