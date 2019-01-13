% revised 03/26/13

function xyt_final = centroids1(file, sz, threshold)

filesi = dir(file); % lists files

cnt = [];
xyt = [];

% particle peak find
for i=1 % loop over each image within a time point
    disp(['n=' num2str(i)]) % displays progress
    im = imread([filesi(i).name]);
    im = (double(im)); % double precision for centroid function
    im = im(:,:,1); % just using one channel of RGB
    im_filter = bpass(im,1,sz); % runs filter bpass and saves a new image

    mask = im_filter > threshold; % sets a threshold

    im_binary = bwlabel(mask); % makes a binary image of the masked image and labels each blob numerically

    os_struct = regionprops(im_binary,im,'WeightedCentroid');
    c = reshape([os_struct.WeightedCentroid],2,length(os_struct))';
    
%     os = []; % applies the mask using the blobs
%     for blob = 1:max(max(im_binary))
%         mask2 = im_binary == blob;
%         os = [os; orientation(im,mask2)]; % finds the subpixel x and y of the center of each blob
%     end

%     c = [os.x;os.y]'; % saving the x and y to matrix c
    c(:,3) = i; % assigns the image number (time) to the last column of c

    cnt = [cnt; c]; % appends c to cnt

end

xyt(:,1) = cnt(:,1); % assigns the x-coordinate to 1st column of xyt
xyt(:,2) = cnt(:,2); % assigns the y-coordinate to 2nd column of xyt
xyt(:,3) = cnt(:,3); % assigns the frame number to 3rd column of xyt

xyt_final = xyt;
