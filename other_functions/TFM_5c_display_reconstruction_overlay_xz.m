disp('Displaying the reconstructed field');
tic
% Px Pz Py Tx Tz Ty are the components of P and T, and PxC PyC PzC TxC TyC
% TzC are the components of P and T with nan values when the CC is under
% the threshold value. Indexes are (i,j,k,t)
if ~plot_force
    UxFit=Ufit(:,:,:,:,1);
    UyFit=Ufit(:,:,:,:,2);
    UzFit=Ufit(:,:,:,:,3);
    UxFit(segmented_cell)=nan;
    UyFit(segmented_cell)=nan;
    UzFit(segmented_cell)=nan;
end
% Similarly for the force:
if plot_force
    FxFit=Ffit(:,:,:,:,1);
    FyFit=Ffit(:,:,:,:,2);
    FzFit=Ffit(:,:,:,:,3);
end
% Extend original image to include one more channel (e.g. blue arrows) and resize: 
img4=permute(imresize(permute(cat(3,img(:,:,:,:,timepoints_to_process),zeros(size(img(:,:,1,:,timepoints_to_process)))),[1,4,3,2,5]),f_pic),[1,4,3,2,5]);
%img4=img4(1:floor(n_px(1)*f_pic),:,:,1:floor(n_layers*f_pic),:);

% Then scan through the layers for which I have a vector field, and plot the vector field in the size of the picture:
clear overlay
t2=1;
for t=timepoints_to_process
    for j=1:ny
        y=yrange(j);
        hold off
        iptsetpref('ImshowBorder','tight')
        %imshow(cat(3,img_ref(:,:,z,t),img_t(:,:,z),zeros(n_px(1),n_px(2))));
        tmp=size(img4);
        imshow(zeros(tmp(1),tmp(4)));
        hold on
        if plot_force
            quiver(Pz(:,j,:,t)*f_pic,Px(:,j,:,t)*f_pic,FzFit(:,j,:,t)*f_F*f_pic,FxFit(:,j,:,t)*f_F*f_pic,0,'color','b');
        else
            quiver(Pz(:,j,:,t)*f_pic,Px(:,j,:,t)*f_pic,UzFit(:,j,:,t)/pxsize_z*f_U*f_pic,UxFit(:,j,:,t)/pxsize_xy*f_U*f_pic,0,'color','b');
        end
        
        tmp = getframe(gcf);
        tmp=tmp.cdata;  %xyc
        tmp=tmp(:,:,3); %arrows only
        if size(tmp,2)>size(img4,4)     % If the image is very narrow in the z direction, margins are created to fit the width of the menus, needs to be removed. 
            margin=floor((size(tmp,2)-size(img4,4))/2);
            tmp=tmp(:,margin+1:margin+size(img4,4));
        end
        img4(:,y,n_colors+1,:,t2)=tmp;
    end
    t2=t2+1;
    close(gcf)
end
img4 = permute(img4,[1 4 3 2 5]);
img4=reshape(img4,size(img4,1),size(img4,2),size(img4,3),size(img4,4)*size(img4,5));
if show_movie_fitted
    if size(img4,3)==3
        implay(img4);
    elseif size(img4,3)==2
        implay(cat(3,img4,zeros(size(img4(:,:,1,:)))));
    end
end
if save_overlay_fitted
    if ~exist('output','dir')
        mkdir('output')
    end
    if single_layer
        if plot_force
            filename='output\overlay_force_reconstruction_xz.tif';
        else
            filename='output\overlay_fitted_displacement_xz.tif';
        end
        imsave3D(img4(:,:,:,layer_saved_xz+n_px(2)*(0:size(timepoints_to_process,2)-1)),filename);
    else
        if plot_force
            filename='output\overlay_force_reconstruction_xz.tif';
        else
            filename='output\overlay_fitted_displacement_xz.tif';
        end
        imsave3D(img4(:,:,:,:),filename);
    end
end
%clear img_4
toc
close(gcf)
