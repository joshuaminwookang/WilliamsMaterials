function [fitresult, gof] = fit_poly22_to_points(r)

[xData, yData, zData] = prepareSurfaceData( r(:,1), r(:,2), r(:,3));

% Set up fittype and options.
ft = fittype( 'poly22' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

if 0
% Plot fit with data.
figure( 'Name', 'poly22 fit' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'poly22 fit', 'z vs. x, y', 'Location', 'NorthEast' );
% Label axes
xlabel( 'x' );
ylabel( 'y' );
zlabel( 'z' );
grid on
view( 1.5, 32.0 );
axis equal
end

