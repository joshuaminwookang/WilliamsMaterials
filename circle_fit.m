% This is a function to fit a circle to a given range of data. Fits for
% radius and center position.


function r0z0R_fit = circle_fit(surface_profile_full,varargin)

if ~isempty(varargin)
    zlimits = varargin{1}(1:2);
    if numel(varargin{1})>2
        xlimits = varargin{1}(3:4);
    else
        xlimits = [Inf -Inf];
    end
else
    zlimits = [Inf -Inf]; 
    xlimits = [Inf -Inf]; %if range is not specified, will include all data in the profile
end

% excise just the sphere profile:
sphere_profile_points = surface_profile_full(:,2) > zlimits(2) & surface_profile_full(:,2) < zlimits(1) ...
                        & surface_profile_full(:,1) > xlimits(2) & surface_profile_full(:,1) < xlimits(1);
sphere_profile_r = surface_profile_full(sphere_profile_points,1);
sphere_profile_z = surface_profile_full(sphere_profile_points,2);

% ...and fit this to a circle...
% ...in r,theta coordinates...

r0z0R_fit = fminsearch(@(x0y0R) circle_fit_error(x0y0R,sphere_profile_r,sphere_profile_z),[mean(sphere_profile_r) mean(sphere_profile_z) max(sphere_profile_r)-mean(sphere_profile_r)]);


if 1 %run_interactively
    figure(72947); hold off
    plot(r0z0R_fit(1)-surface_profile_full(:,1), surface_profile_full(:,2),'.')
    hold all
    plot(r0z0R_fit(1)-sphere_profile_r,sphere_profile_z,'.')
    viscircles([0 r0z0R_fit(2)],r0z0R_fit(3),'EdgeColor','k');
    
grid on
axis equal
%ylim([-20 40]); xlim([-40 40])
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
xlabel('R (µm)')
ylabel('Surface displacement (µm)')
pause(0.5)
end
