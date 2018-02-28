disp('Display the overlay')
tic

% If only rendering one layer, pick up the closest layer with known data:
if single_layer
    if layer_saved<1
        layer_saved_tmp=layer_saved*max(zrange);
    end
    tmp=zrange(2)-zrange(1);
    k=find(zrange>=layer_saved_tmp-tmp/2,1);
    layer_saved_tmp=zrange(k);
end

% Then extend original image to include one more channel (e.g. blue arrows) and resize:
if single_layer
    img4=imresize(cat(3,img(:,:,:,layer_saved_tmp,timepoints_to_process),zeros(size(img(:,:,1,layer_saved_tmp,timepoints_to_process)))),f_pic);
    img4=img4(1:floor(n_px(1)*f_pic),1:floor(n_px(2)*f_pic),:,1,:);
else
    img4=imresize(cat(3,img(:,:,:,:,timepoints_to_process),zeros(size(img(:,:,1,:,timepoints_to_process)))),f_pic);
    img4=img4(1:floor(n_px(1)*f_pic),1:floor(n_px(2)*f_pic),:,:,:);
end

% Then scan through the layers for which I have a vector field, and plot the vector field in the size of the picture:
clear overlay
t2=1;
figure
for t=timepoints_to_process
    if ~single_layer
        for k=1:nz
            z=zrange(k);
            hold off
            iptsetpref('ImshowBorder','tight')
            %imshow(cat(3,img_ref(:,:,z,t),img_t(:,:,z),zeros(n_px(1),n_px(2))));
            imshow(zeros(floor(n_px(1)*f_pic),floor(n_px(2)*f_pic)));
            hold on
            quiver(PyC(:,:,k,t)*f_pic,PxC(:,:,k,t)*f_pic,UyC(:,:,k,t)/pxsize_xy*f_U*f_pic,UxC(:,:,k,t)/pxsize_xy*f_U*f_pic,0,'color','b');
            tmp = getframe(gcf);
            tmp=tmp.cdata;  %xyc
            tmp=tmp(:,:,3); %arrows only
            img4(:,:,n_colors+1,z,t2)=tmp;
        end
        t2=t2+1;
    else
        hold off
        iptsetpref('ImshowBorder','tight')
        %imshow(cat(3,img_ref(:,:,z,t),img_t(:,:,z),zeros(n_px(1),n_px(2))));
        imshow(zeros(floor(n_px(1)*f_pic),floor(n_px(2)*f_pic)));
        hold on
        quiver(PyC(:,:,k,t)*f_pic,PxC(:,:,k,t)*f_pic,UyC(:,:,k,t)/pxsize_xy*f_U*f_pic,UxC(:,:,k,t)/pxsize_xy*f_U*f_pic,0,'color','b');
        tmp = getframe(gcf);
        tmp=tmp.cdata;  %xyc
        tmp=tmp(:,:,3); %arrows only
        img4(:,:,n_colors+1,1,t2)=tmp;
        t2=t2+1;
    end
end
close(gcf);
img4=reshape(img4,size(img4,1),size(img4,2),size(img4,3),size(img4,4)*size(img4,5));

if show_movie_experimental
    if size(img4,3)==3
        implay(img4);
    elseif size(img4,3)==2
        implay(cat(3,img4,zeros(size(img4(:,:,1,:)))));
    end
end

if save_overlay_experimental
    disp('Saving the overlay');
    if ~exist('output','dir')
        mkdir('output')
    end
    
    if single_layer
        filename='output\overlay_experimental_displacement_single_layer.tif';
        img4=img4(:,:,:,:);
        imsave3D(img4(:,:,:,:),filename);
    else
        filename='output\overlay_experimental_displacement.tif';
        imsave3D(img4,filename);
    end
    disp(filename);
end

clear img_4
close(gcf)
