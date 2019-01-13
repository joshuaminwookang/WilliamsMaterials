%show_slice(slices, r, imagestack)
%this is a function to visually check the results of 3D particle locating. 
%the variable "slices" should be a vector list of the slices you want to look at;
%for example, if you wanted to look at slices 1 through 30 of a stack, you should have 
%slices = [1:30]; 
%if you wanted just to check every other one, you could do like 
%slices = [1:2:30]
%
%r is the matrix of particle coordinates with the x-, y-, and z-coordinates
%in the first 3 columns (standard output of Kilfoil's tracking). Should be
%in PIXELS (so the numbers match up with the imagestack size)
%
%imagestack is the 3-dimensional array; standard output of Kilfoil's
%bandpass3dMB program OR the raw image stack (often called XYZ for us)
%
%if you want to check, for example, the X-Z tracking, you should just
%permute things around when you call it, e.g. :
% show_slice([1:250], r(:,[3 2 1]), permute(imagestack,[3 2 1])) %<-- this yields an X-Z stack

function show_slice(slices, r, imagestack)
sorted = sortrows(r,3);
figure
colormap jet


k=1;
while k<=numel(slices)
    display(['slice: ' num2str(k)])
    dz = 4.2; %<-- this is how many steps +/- you consider to still be "in the layer", roughly your particle radius in pixels -- depending on your pixel size, you might want to adjust this
    inlayer = (sorted(:,3)>slices(k)-dz)&(sorted(:,3)<slices(k)+dz);
    %inlayer = (sorted(:,3)>0.1497*slices(k)-dz)&(sorted(:,3)<0.1497*slices(k)+dz);
    nearlayer_x = (sorted(inlayer,1));
    nearlayer_y = (sorted(inlayer,2));
    
    
    imagesc(imagestack(:,:,slices(k))); hold on; plot(nearlayer_y,nearlayer_x,'ko','MarkerSize',8,'MarkerFaceColor','w','LineWidth',2); hold off
    %axis image; %axis([100 170 100 170])
    axis xy; axis equal; axis tight
    title(k)
    
    %%commented out below are some examples of manually putting in the
    %%scaling so that your axes are in microns (the first one is for our
    %%silica particles)
    
    %imagesc([0:512]*(91.18/512), [0:512]*(91.18/512),imagestack(:,:,slices(k))); hold on; plot(nearlayer_y*91.18/512,nearlayer_x*91.18/512,'k+'); hold off
    %imagesc([1:256]*(0.2062), [1:256]*(0.2062),imagestack(:,:,slices(k))); hold on; plot(nearlayer_y,nearlayer_x,'k+'); hold off
    
   
    

    %axis xy %<--- makes the X-Z version put z=0 at the bottom of the picture
    
    i=input('Press Enter to check the next slice on the list. Type ''b'', then Enter, to go back a slice.','s');
    
    if strcmp(i,'b')
        k=k-1;
    else
        k=k+1;
    end
    
end