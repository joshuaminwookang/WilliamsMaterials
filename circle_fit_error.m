%the error in a fit to a circle ... need to fit center position (x0,y0) and
%radius, R
function [calc_error, residuals] = circle_fit_error(x0y0R,x,y)

x0 = x0y0R(1);
y0 = x0y0R(2);
R = x0y0R(3);

%convert the (x,y) data to radial distance from (x0, y0)
R_data = sqrt(sum([(x-x0).^2 (y-y0).^2],2));

%what's the error with respect to the current best-guess for the radius?
calc_error = sqrt(sum((R_data - R).^2));
residuals = R_data - R;

end
    