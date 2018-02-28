% If recomputing displacements, the following force reconstruction
% variables will have to be reinitialized if already computed:
clear Ufit Ffit cell_surface_found

% box_size2 can be simply divided by 2 to have the number of
% values in the box on each side of the center point. 
box_size_xy2=box_size_xy-1;
box_size_z2=box_size_z-1;

% Separating layers and timepoints, and zero-padding:
border_xy=box_size_xy2/2+max_distance_correlation_xy+2;
border_z=box_size_z2/2+max_distance_correlation_z+2;

img_padded=zeros([n_px+2*border_xy,n_layers+2*border_z,n_timepoints]);
for t=1:n_timepoints
    img_padded(:,:,:,t)=padarray(squeeze(img(:,:,tracking_channel,:,t)),[border_xy border_xy border_z]);
end

% Redefine the dimensions taking into account the pads: 
sizeimg_padded=size(img_padded);
n_px2=sizeimg_padded(1:2);
n_layers2=sizeimg_padded(3);

% Size of the reference box: should fit the correlation distance + the size of the moving box. 
hw=max_distance_correlation_xy+box_size_xy2/2;
hh=max_distance_correlation_z+box_size_z2/2;

disp('Cross-correlation');
tic
% Range used in the padded images:
xrange2=border_xy+box_size_xy2/2+1:step_size_xy:n_px2(1)-border_xy-box_size_xy2/2-1;   % The first value in the range corresponds to the first box fitting in the actual picture out of the zero pad. 
yrange2=border_xy+box_size_xy2/2+1:step_size_xy:n_px2(2)-border_xy-box_size_xy2/2-1;
zrange2=border_z+box_size_z2/2+1:step_size_z:n_layers2-border_z-box_size_z2/2-1;

% Equivalent ranges without the pad (used for the display on original picture): 
xrange=box_size_xy2/2+1:step_size_xy:n_px(1)-box_size_xy2/2-1;   % The first value in the range corresponds to the first box fitting in the actual picture out of the zero pad. 
yrange=box_size_xy2/2+1:step_size_xy:n_px(2)-box_size_xy2/2-1;
zrange=box_size_z2/2+1:step_size_z:n_layers-box_size_z2/2-1;

nx=size(xrange2,2);
ny=size(yrange2,2);
nz=size(zrange2,2);
disp('    nx   ny   nz');
disp([nx ny nz]);
disp('     z    nz    t     n_timepoints');
Uexp=zeros(nx,ny,nz,n_timepoints,3); % Translation vectors, in pixels, for each position considered. 
P=zeros(nx,ny,nz,n_timepoints,3); % Position vectors, in pixels, where Ts are computed. 
C=zeros(nx,ny,nz,n_timepoints);   % Correlation values at the optimum. 

fftw('planner', 'measure');

ti=1;
if reference_frame==-1
    ti=2;
    img_ref=img_padded(:,:,:,1);
else
    img_ref=img_padded(:,:,:,reference_frame);
end

for t=ti:n_timepoints
    if t~=reference_frame 
        if t>2 && reference_frame==-1
            img_ref=img_t;
        end
        img_t=img_padded(:,:,:,t);

        parfor k=1:nz   % xyz is the position of the center of the box that is carried around. Leave a half box size out on the side (so whole box until the center) so the box can be translated without going out, +1 pixel so even with the subpixel step it doesnt go out. 
%        for k=1:nz   % xyz is the position of the center of the box that is carried around. Leave a half box size out on the side (so whole box until the center) so the box can be translated without going out, +1 pixel so even with the subpixel step it doesnt go out. 
            disp([k nz t n_timepoints]);
            z=zrange2(k);
            for j=1:ny
                y=yrange2(j);
                for i=1:nx
                    x=xrange2(i);
                    box_ref=double(img_ref(x-hw:x+hw,y-hw:y+hw,z-hh:z+hh));
                    box=img_t(x-box_size_xy2/2:x+box_size_xy2/2,y-box_size_xy2/2:y+box_size_xy2/2,z-box_size_z2/2:z+box_size_z2/2);
                    TxyzCC=xcorr3subpxnormQuadInt(box,box_ref);
                    Txyz=TxyzCC(1:3);                                           % T is the image translation in pixels
                    CC=TxyzCC(4);
                    P(i,j,k,t,:)=[x y z]-[border_xy-1 border_xy-1 border_z-1];
                    Uexp(i,j,k,t,:)=[Txyz(1:2)*pxsize_xy, Txyz(3)*pxsize_z];    % U is the physical displacement in um
                    C(i,j,k,t)=CC;
                end
            end
        end
    end
end
clear img_padded img_t img_ref
