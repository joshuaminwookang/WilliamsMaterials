% This is a function that returns a logical vector of length N (number of particles) where the ones 
% correspond to the particles that are included in the selected region of
% interest. This allows the user to designate a desired geometric set of
% particles while keeping the remainder of the particles' positions in
% memory for use in calculations.
%
% function select_vector = select_region(XYZ, shape, varargin)
%
% Function expects input 'XYZ' to be a 2-d array of length N with the XYZ particle positions
% in columns 1:3.
%
% 'shape' must be a string designating the type of shape selection.
% 
% select_vector = select_region(XYZ, 'sphere', radius_from_center); 
% will select a spherical region beginning in the geometric center of XYZ 
%
% select_vector = select_region(XYZ, 'sphere', radius_from_center, center_point); 
% will select a spherical region beginning at the coordinates center_point
%
% select_vector = select_region(XYZ, 'relative', xmargin, ymargin, zmargin);
% will select a rectangular region smaller than the extremes of the XYZ
% sample by the designated margins
%
% select_vector = select_region(XYZ, 'relative', margin);
% will select a rectangular region smaller than the extremes of the XYZ
% sample by the same margin in the x-,y-, and z-directions
%
% select_vector = select_region(XYZ, 'box', minima, maxima);
% will select a rectangular box region with diagonally opposite corners at
% minima = [xmin ymin zmin] and maxima = [xmax ymax zmax]
%
% Kate Jensen - May 26, 2009; finished July 31, 2009

function select_vector = select_region(XYZ, shape, varargin)

%start with the box shape, since I expect to use that the most

if strcmp(shape,'relative')
    %this lets the user specify relative bounds
    
    xmargin = varargin{1};
    
    if length(varargin) == 1 %if user only inputs one number, use same value for x,y,and z-margins
        ymargin = xmargin;
        zmargin = xmargin;
    else
        ymargin = varargin{2};
        zmargin = varargin{3};
    end
        
    maxima = max(XYZ);
    minima = min(XYZ);
    
    select_vector = ...
        (XYZ(:,1)<=(maxima(1)-xmargin) & XYZ(:,1)>=(minima(1)+xmargin)) &...
        (XYZ(:,2)<=(maxima(2)-ymargin) & XYZ(:,2)>=(minima(2)+ymargin)) &...
        (XYZ(:,3)<=(maxima(3)-zmargin) & XYZ(:,3)>=(minima(3)+zmargin));
end

if strcmp(shape,'box')
    %this lets the user specify particular bounds
    
    minima = varargin{1};
    maxima = varargin{2};
    
    select_vector = ...
        (XYZ(:,1)<=maxima(1) & XYZ(:,1)>=minima(1)) &...
        (XYZ(:,2)<=maxima(2) & XYZ(:,2)>=minima(2)) &...
        (XYZ(:,3)<=maxima(3) & XYZ(:,3)>=minima(3));
    
    
end


if strcmp(shape,'sphere')
    radius_from_center = varargin{1};
    %let user specify a starting center point; otherwise use the sample center
    if length(varargin) == 1    
        maxima = max(XYZ);
        minima = min(XYZ);
        center_point = minima+(maxima-minima)/2;
    else
        center_point = varargin{2};
    end

    
    center_point = ones(length(XYZ),1)*center_point;
    dist_from_center = sqrt(sum((XYZ-center_point).^2,2));
    
    select_vector = dist_from_center <= radius_from_center;    
end

