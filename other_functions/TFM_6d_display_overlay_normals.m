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
img4=imresize(cat(3,img(:,:,:,:,timepoints_to_process),zeros(size(img(:,:,1,:,timepoints_to_process)))),f_pic);
img4=img4(1:floor(n_px(1)*f_pic),1:floor(n_px(2)*f_pic),:,:,:);

% Then scan through the layers for which I have a vector field, and plot the vector field in the size of the picture:
clear overlay
t2=1;
zspacing=P_fraction*(zrange(2)-zrange(1))/2;
figure;
for t=timepoints_to_process
    for k=1:nz
        z=zrange(k);
        hold off
        iptsetpref('ImshowBorder','tight')
        %imshow(cat(3,img_ref(:,:,z,t),img_t(:,:,z),zeros(n_px(1),n_px(2))));
        imshow(zeros(floor(n_px(1)*f_pic),floor(n_px(2)*f_pic)));
        hold on
        
        % Find the range associated with each z position:
        if k==1
            associated_range=[z z+zspacing]*pxsize_z;
        elseif k==nz
            associated_range=[z-zspacing z]*pxsize_z;
        else
            associated_range=[z-zspacing z+zspacing]*pxsize_z;
        end
        
        % Pick up the normals vectors in this range:
        indices=find(logical(and(positions{t}(:,3)>associated_range(1), positions{t}(:,3)<associated_range(2))));
        positions_tmp=positions{t}(indices,:);
        normals_vectors_tmp=normals{t}(indices,:);
        
        % Plot the vector field:
        quiver(positions_tmp(:,2)/pxsize_xy*f_pic,positions_tmp(:,1)/pxsize_xy*f_pic,normals_vectors_tmp(:,2)/pxsize_xy*15*f_pic,normals_vectors_tmp(:,1)/pxsize_xy*15*f_pic,0,'color','b');
        
        % Uncomment to place a star at the contraction center (forces pointed towards it):
        %plot(O(1,t),O(2,t),'*');
        tmp = getframe(gcf);
        tmp=tmp.cdata;  %xyc
        tmp=tmp(:,:,3); %arrows only
        img4(:,:,n_colors+1,z,t2)=tmp;
    end
    t2=t2+1;
end
close(gcf);
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
    filename='output\overlay_normals_stack.tif';
    imsave3D(img4(:,:,:,:),filename);
end
clear img_4
toc