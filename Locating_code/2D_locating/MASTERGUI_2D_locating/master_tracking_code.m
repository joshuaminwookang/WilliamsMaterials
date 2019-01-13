% Code that takes an image and a structure from
% master_tracking_gui, and then applies the tracking parameters from the
% structure to the image
% Written by Rob Dec 2016


function pks=master_tracking_code(imagefile,tracking_struct)
pks=[];
if isa(imagefile,'char')==1
    im=imread(imagefile);
else
    im=imagefile;
end

% Do the bpass
if isempty(tracking_struct.bpass_params)
    error('No bpass parameters supplied')
elseif length(tracking_struct.bpass_params)~=3
    error('Wrong number of bpass parameters supplied (should be 3)')
else
    im_bpass=bpass(double(im),tracking_struct.bpass_params(1),tracking_struct.bpass_params(2),tracking_struct.bpass_params(3));
end

if isempty(tracking_struct.roughfind_method)
    error('No rough find method supplied')
elseif strcmp(tracking_struct.roughfind_method,'pkfnd')
    % PUT PKFND STUFF HERE
    if isempty(tracking_struct.roughfind_params)
        error('No pkfnd parameters supplied')
    elseif length(tracking_struct.roughfind_params)~=2
        error('Wrong number of pkfnd parameters supplied (should be 2)')
    else
        p=pkfnd(im_bpass,tracking_struct.roughfind_params(1),tracking_struct.roughfind_params(2));
    end
elseif strcmp(tracking_struct.roughfind_method,'binarize')
    % PUT STUFF HERE ABOUT BINARIZING
    if isempty(tracking_struct.roughfind_params)
        error('No binarize parameters supplied')
    elseif length(tracking_struct.roughfind_params)~=3
        error('Wrong number of binarize parameters supplied (should be 3)')
    else
        BW=binarize(im_bpass,tracking_struct.roughfind_params(1),tracking_struct.roughfind_params(2),tracking_struct.roughfind_params(3));
        regions = bwlabel(BW,8);
        %props = regionprops(regions, 'Centroid', 'Area'); % find features in BW (Particles)
        props = regionprops(regions,im,'WeightedCentroid','Area'); % find features in BW (Particles)

        X=[]; Y=[]; Area=[];
        for n = 1:1:size(props)    % go through features in BW
            % Particles:
                X = [X; props(n).WeightedCentroid(1)];
                Y = [Y; props(n).WeightedCentroid(2)];
                Area = [Area; props(n).Area];   
        end

        pks=cat(2,X,Y,Area);
    end
else
    error('Unrecognised roughfind_method')
end

if isempty(pks)
   if isempty(tracking_struct.refine_method)
       error('No refine_method supplied')
   elseif strcmp(tracking_struct.refine_method,'cntrd')
       if length(tracking_struct.refine_params)~=1
           error('Unexpected number of refine_parameters (should be 1)')
       else
           pks=cntrd(im_bpass,p,tracking_struct.refine_params);
       end
   elseif strcmp(tracking_struct.refine_method,'radialcenter')
       if length(tracking_struct.refine_params)~=1
           error('Unexpected number of refine_parameters (should be 1)')
       else
               particle_sz=tracking_struct.refine_params;
                pks=[];
                %zeros(length(p),3);
                for i=1:length(p)
                    xcentre=round(p(i,1)); % Centre of the particle to pixel level resolution
                    ycentre=round(p(i,2));
                    if ycentre-(particle_sz-1)/2>0 && ycentre+(particle_sz-1)/2<=length(im(:,1)) && xcentre-(particle_sz-1)/2 >0 && xcentre+(particle_sz-1)/2 <= length(im(1,:))
                        % Here we take an area of the image that is particle_sz by
                        % particle_sz in dimensions around the position of each
                        % particle. Then we find the centre of radial symmetry in this
                        % test area.
                        testim=im(ycentre-(particle_sz-1)/2:ycentre+(particle_sz-1)/2,xcentre-(particle_sz-1)/2:xcentre+(particle_sz-1)/2);
                        [c,d,f]=radialcenter(double(testim));
                        pks=[pks ; xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
                        %p_rs(i,:)=[xcentre-(particle_sz+1)/2+c , ycentre-(particle_sz+1)/2+d , f];
                    end
                end         
       end
   else
       error('Unrecognised refine_method')
   end
end 

    
end
