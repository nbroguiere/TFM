disp('Displaying the reconstructed field');
tic

if single_layer
    % Px Pz Py Tx Tz Ty are the components of P and T, and PxC PyC PzC TxC TyC
    % TzC are the components of P and T with nan values when the CC is under
    % the threshold value. Indexes are (i,j,k,t)
    % Pick up the closest layer with known data:
    if layer_saved<1
        layer_saved_tmp=layer_saved*max(zrange);
    end
    tmp=zrange(2)-zrange(1);
    k=find(zrange>=layer_saved_tmp-tmp/2,1);
    layer_saved_tmp=zrange(k);

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
    % Extend original image to include one more channel (e.g. blue arrows), pick up the layer to show, and resize: 
    img4=squeeze(cat(3,img(:,:,:,layer_saved_tmp,timepoints_to_process),zeros(size(img(:,:,1,1,timepoints_to_process)))));
    img4=imresize(img4,f_pic);
    img4=img4(1:floor(n_px(1)*f_pic),1:floor(n_px(2)*f_pic),:,:);

    % Then scan through the layers for which I have a vector field, and plot the vector field in the size of the picture:
    clear overlay
    i=1;
    for t=timepoints_to_process
            z=layer_saved;
            hold off
            iptsetpref('ImshowBorder','tight')
            imshow(zeros(floor(n_px(1)*f_pic),floor(n_px(2)*f_pic)));
            hold on
            if plot_force
                quiver(Py(:,:,k,t)*f_pic,Px(:,:,k,t)*f_pic,FyFit(:,:,k,t)*f_F*f_pic,FxFit(:,:,k,t)*f_F*f_pic,0,'color','b');
            else
                quiver(Py(:,:,k,t)*f_pic,Px(:,:,k,t)*f_pic,UyFit(:,:,k,t)/pxsize_xy*f_U*f_pic,UxFit(:,:,k,t)/pxsize_xy*f_U*f_pic,0,'color','b');
            end
            % Uncoment to place a star at the contraction center
            %plot(O(1,t),O(2,t),'*');
            tmp = getframe(gcf);
            tmp=tmp.cdata;  %xyc
            tmp=tmp(:,:,3); %arrows only

            img4(:,:,n_colors+1,i)=tmp;
            i=i+1;
    end
    close(gcf)

    if show_movie_fitted
        if length(timepoints_to_process)==1
            imshow(img4)
        else
            if size(img4,3)==3
                implay(img4);
            elseif size(img4,3)==2
                implay(cat(3,img4,zeros(size(img4(:,:,1,:)))));
            end
        end
    end
    if save_overlay_fitted
        if ~exist('output','dir')
            mkdir('output')
        end
        if plot_force
            filename='output\overlay_force_reconstruction_single_layer.tif';
        else
            filename='output\overlay_fitted_displacement_single_layer.tif';
        end
        imsave3D(img4(:,:,:,:),filename);
    end
    clear img_4
else
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
    for t=timepoints_to_process
        for k=1:nz
            z=zrange(k);
            hold off
            iptsetpref('ImshowBorder','tight')
            %imshow(cat(3,img_ref(:,:,z,t),img_t(:,:,z),zeros(n_px(1),n_px(2))));
            imshow(zeros(floor(n_px(1)*f_pic),floor(n_px(2)*f_pic)));
            hold on
            if plot_force
                quiver(Py(:,:,k,t)*f_pic,Px(:,:,k,t)*f_pic,FyFit(:,:,k,t)*f_F*f_pic,FxFit(:,:,k,t)*f_F*f_pic,0,'color','b');
            else
                quiver(Py(:,:,k,t)*f_pic,Px(:,:,k,t)*f_pic,UyFit(:,:,k,t)/pxsize_xy*f_U*f_pic,UxFit(:,:,k,t)/pxsize_xy*f_U*f_pic,0,'color','b');
            end
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
        if plot_force
            filename='output\overlay_force_reconstruction_stack.tif';
        else
            filename='output\overlay_fitted_displacement_stack.tif';
        end
        imsave3D(img4(:,:,:,:),filename);
    end
    clear img_4
end
toc