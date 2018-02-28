% Segment the cell and set cross-correlation to 0 for the elements which are on top of the cell:
tic
clear img_threshold_smooth      % If redoing the segmentation, this will have to be recalculated as well in the display function. 

% Find nx ny nz in the same way as used in the find displacement function
% (dimensions of the mesh on which the mechanical model is solved, needed
% for segmentation of the cell surface):
box_size_xy2=box_size_xy-1;
box_size_z2=box_size_z-1;
border_xy=box_size_xy2/2+max_distance_correlation_xy+2;
border_z=box_size_z2/2+max_distance_correlation_z+2;
xrange=box_size_xy2/2+1:step_size_xy:n_px(1)-box_size_xy2/2-1;   % The first value in the range corresponds to the first box fitting in the actual picture out of the zero pad. 
yrange=box_size_xy2/2+1:step_size_xy:n_px(2)-box_size_xy2/2-1;
zrange=box_size_z2/2+1:step_size_z:n_layers-box_size_z2/2-1;
nx=length(xrange);
ny=length(yrange);
nz=length(zrange);

disp('Cell smoothing and thresholding')
img_threshold=zeros(n_px(1),n_px(2),n_layers,n_timepoints,'logical');   % x y z t
segmented_cell=zeros([nx,ny,nz,n_timepoints],'logical');
if segment_reference_frame
    timepoints_to_segment=1:n_timepoints;
else
    timepoints_to_segment=[1:reference_frame-1 reference_frame+1:n_timepoints];
end
for t=timepoints_to_segment
    % Pick up the stack for the timepoint under consideration:
    img_threshold_tmp=double(squeeze(img(:,:,channel_cell,:,t)))/255;  
    
    % Smooth it:
    box_range=-floor(3*sigma/pxsize_xy):floor(3*sigma/pxsize_xy);
    filter=normpdf(box_range,0,sigma/pxsize_xy);
    img_threshold_tmp=imfilter(img_threshold_tmp,filter);
    img_threshold_tmp=imfilter(img_threshold_tmp,filter');
    box_range=-floor(3*sigma/pxsize_z):floor(3*sigma/pxsize_z);
    filter=zeros(1,1,length(box_range));
    filter(:)=normpdf(box_range,0,sigma/pxsize_z);
    img_threshold_tmp=imfilter(img_threshold_tmp,filter);
    
    % Threshold it:
    img_threshold(:,:,:,t)=img_threshold_tmp>(threshold*max(img_threshold_tmp(:)));
end

disp('Morphological cleaning operations')
for t=timepoints_to_segment
    % Closing
    if floor(2*closing_radius/pxsize_z)>0
        if fast_segmentation
            se=strel('cuboid',[floor(2*closing_radius/pxsize_xy),floor(2*closing_radius/pxsize_xy),floor(2*closing_radius/pxsize_z)]);
        else
            se=strel('sphere',floor(closing_radius/pxsize_xy));
        end
        img_threshold(:,:,:,t)=imclose(img_threshold(:,:,:,t),se);
    end
    % Filling holes
    if fill_holes
        for z=1:n_layers
            img_threshold(:,:,z,t)=imfill(img_threshold(:,:,z,t),'holes');
        end
    end
    % Keep only the biggest object (e.g. to remove debris but keep the main cell). 
    if keep_only_biggest
        if tolerance_keep_biggest>0
            if fast_segmentation
                se=strel('cuboid',[floor(2*tolerance_keep_biggest/pxsize_xy),floor(2*tolerance_keep_biggest/pxsize_xy),floor(2*tolerance_keep_biggest/pxsize_z)]);
            else
                se=strel('sphere',floor(tolerance_keep_biggest/pxsize_xy));
            end
            img_threshold_tmp=imdilate(img_threshold(:,:,:,t),se);
        else
            img_threshold_tmp=img_threshold(:,:,:,t);
        end
        CC=bwconncomp(img_threshold_tmp);
        stats=regionprops(CC,'area');
        [value,pos]=max([stats.Area]);
        indices=CC.PixelIdxList{pos};
        mask=zeros(size(img_threshold_tmp),'logical');
        mask(indices)=1;
        img_threshold(:,:,:,t)=img_threshold(:,:,:,t).*mask;
        clear CC stats indices mask
    end
    % Mark the position of the segmented cell in the mechanics mesh:
    for i=1:nx
        for j=1:ny
            for k=1:nz
                if img_threshold(xrange(i),yrange(j),zrange(k))
                    segmented_cell(i,j,k,t)=1;
                end
            end
        end
    end
end

clear img_threshold_tmp filter
toc

% Find out which positions are located on the cell surface and inside the cell. 
disp('finding cell surface')
tic
cell_position=zeros([nx,ny,nz,n_timepoints],'logical');    % indices: i j k t
cell_surface=zeros([nx,ny,nz,n_timepoints],'logical');     % indices: i j k t
sxyp=floor(step_size_xy/2);
sxym=floor((step_size_xy-1)/2);
szp=floor(step_size_z/2);
szm=floor((step_size_z-1)/2);
numelbox=step_size_xy.^2*step_size_z;
% Definition of the neighborhood used to find the surface (e.g. 6 neighbourhood)
neighborhood=zeros([3,3,3],'logical');
neighborhood(1:3,2,2)=1;
neighborhood(2,1:3,2)=1;
neighborhood(2,2,1:3)=1;
for t=timepoints_to_segment
    for k=1:nz
        z=zrange(k);
        zmin=max(z-szm,1);
        zmax=min(z+szp,n_layers);
        for j=1:ny
            y=yrange(j);
            ymin=max(y-sxym,1);
            ymax=min(y+sxyp,n_px(2));
            for i=1:nx
                x=xrange(i);
                xmin=max(x-sxym,1);
                xmax=min(x+sxyp,n_px(1));
                tmp=img_threshold(xmin:xmax,ymin:ymax,zmin:zmax,t);
                tmp2=sum(double(tmp(:)));
                if tmp2>0
                    cell_position(i,j,k,t)=1;
                end
            end
        end
    end
    for i=1:nx
        for j=1:ny
            for k=1:nz
                if cell_position(i,j,k,t)
                    tmp=cell_position(max(i-1,1):min(i+1,nx),max(j-1,1):min(j+1,ny),max(k-1,1):min(k+1,nz),t);
                    i1=1;i2=3;j1=1;j2=3;k1=1;k2=3;
                    if i==1
                        i1=2;
                    end
                    if i==nx
                        i2=2;
                    end
                    if j==1
                        j1=2;
                    end
                    if j==ny
                        j2=2;
                    end
                    if k==1
                        k1=2;
                    end
                    if k==nz
                        k2=2;
                    end
                    tmp2=neighborhood(i1:i2,j1:j2,k1:k2);
                    tmp3=tmp.*tmp2;
                    if sum(tmp3(:))<sum(tmp2(:))
                        cell_surface(i,j,k,t)=true;
                    end
                end
            end
        end
    end
end
tmp=zeros(nx-2,ny-2,nz-2);
borders=logical(padarray(tmp,[1,1,1],1));
cell_position=xor(cell_position,cell_surface);
toc
