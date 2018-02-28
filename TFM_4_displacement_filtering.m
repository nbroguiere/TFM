%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the experimental displacement Uexp, of components Ux Uy Uz, 
% and computes in return the corrected experimental displacement UexpC, 
% of components UxC UyC UzC. Add NaN values when CC is under the threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fill in the position information for the frame which has no measured
% vectors in it:
if reference_frame>1
    P(:,:,:,reference_frame,:)=P(:,:,:,1,:);
elseif reference_frame==1
    P(:,:,:,reference_frame,:)=P(:,:,:,2,:);
end

% Px Pz Py Tx Tz Ty are the components of P and T
Px=P(:,:,:,:,1);
Py=P(:,:,:,:,2);
Pz=P(:,:,:,:,3);
Ux=Uexp(:,:,:,:,1);
Uy=Uexp(:,:,:,:,2);
Uz=Uexp(:,:,:,:,3);

% Set cross-correlation to 0 for the elements which are on top of the cell:
tic
disp('Computing the rejection area around the cell');
for t=timepoints_to_segment
    % Set the cross-correlation to 0 on the cell, also keep the position of 
    % the cell for later use in the fitted field overlay:
    if rejection_radius>0
        if fast_segmentation
            se=strel('cuboid',[floor(2*rejection_radius/pxsize_xy),floor(2*rejection_radius/pxsize_xy),floor(2*rejection_radius/pxsize_z)]);
        else
            se=strel('sphere',floor(rejection_radius/pxsize_xy));
        end
        img_threshold_tmp=imdilate(img_threshold(:,:,:,t),se);
    else
        img_threshold_tmp=img_threshold(:,:,:,t);
    end
    for i=1:nx
        for j=1:ny
            for k=1:nz
                if img_threshold_tmp(xrange(i),yrange(j),zrange(k))
                    C(i,j,k,t)=0;
                end
            end
        end
    end
end
clear img_threshold_tmp filter
toc

disp('Filtering of the experimental displacement field U')
tic
% PxC PyC PzC UxC UyC UzC are the corrected components of P and U with nan values 
% when the CC is under the threshold value. Indices are (i,j,k,t)
if reference_frame>0
    C(:,:,:,reference_frame)=ones(size(C(:,:,:,reference_frame)));
end
tmp=find(C<autocorrelation_threshold);
UxC=Ux;
UyC=Uy;
UzC=Uz;
UxC(tmp)=nan(size(tmp));
UyC(tmp)=nan(size(tmp));
UzC(tmp)=nan(size(tmp));
PxC=Px;
PyC=Py;
PzC=Pz;
PxC(tmp)=nan(size(tmp));
PyC(tmp)=nan(size(tmp));
PzC(tmp)=nan(size(tmp));

% Subtract the median in each layer:
if subtract_median_layer
    for t=1:n_timepoints
        for k=1:nz
            UxCtmp=UxC(:,:,k,t);
            UyCtmp=UyC(:,:,k,t);
            UzCtmp=UzC(:,:,k,t);
            midx=median(UxCtmp(~isnan(UxCtmp)));
            midy=median(UyCtmp(~isnan(UxCtmp)));  % put again a find before isnan if it crashes. TCEP
            midz=median(UzCtmp(~isnan(UxCtmp)));
            UxC(:,:,k,t)=UxCtmp-midx;
            UyC(:,:,k,t)=UyCtmp-midy;
            UzC(:,:,k,t)=UzCtmp-midz;
        end
    end
end
clear UxCtmp UyCtmp UzCtmp

% Compute the integrals of the displacement
if compute_integrals
    UxC=cumsum(UxC,4,'omitnan');
    UyC=cumsum(UyC,4,'omitnan');
    UzC=cumsum(UzC,4,'omitnan');
end

if subtract_median_position
    UxC=UxC-repmat(median(UxC,4,'omitnan'),[1 1 1 n_timepoints]);
    UyC=UyC-repmat(median(UyC,4,'omitnan'),[1 1 1 n_timepoints]);
    UzC=UzC-repmat(median(UzC,4,'omitnan'),[1 1 1 n_timepoints]);
end
toc

UexpC=cat(5,UxC,UyC,UzC);
