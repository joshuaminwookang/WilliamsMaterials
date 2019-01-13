% function to do a radial collapse of data manually, until the user is
% satisfied with the "goodness" of the collapse
%
% Kate Jensen - functionalized December 10, 2014

function [center,rz] = radial_collapse_manual(r)

dx = 0.1; 

%classy loops to get the center positions...
still_working = 1;

figure(86038); hold off
scatter3(r(:,1),r(:,2),r(:,3),30,r(:,3))
set(gca,'LineWidth',1,'FontSize',20,'FontWeight','bold'); box on; grid on
axis equal
view(0,90)

center = input('starting approx center: ');

direction = 'x';

while still_working
    
    rz = [sqrt((r(:,1)-center(1)).^2 + (r(:,2)-center(2)).^2) r(:,3)];
    plot(rz(:,1),rz(:,2),'.');  xlim([0 80])
    title(center)
    
    s = input(['Press Enter to increase by ' num2str(dx) ' µm in ' direction ', ''d'' to decrease, specify axis (x/y), or ''f'' for "finished": '],'s');
    
    if strcmp(s,'f')
        %display('all done!')
        still_working = 0;
    elseif strcmp(s,'x') || strcmp(s,'y')
        direction = s;
    elseif strcmp(s,'d')
        if strcmp(direction,'x')
            center(1) = center(1) - dx;
        else
            center(2) = center(2) - dx;
        end
    else
        if strcmp(direction,'x')
            center(1) = center(1) + dx;
        else
            center(2) = center(2) + dx;
        end
    end
    
end
