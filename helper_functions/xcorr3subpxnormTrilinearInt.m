function [shiftxyz_CC] = xcorr3subpxnormTrilinearInt(box,ref,subpxres)
% shiftxyz_CC = xcorr3subpxnormTrilinearInt(box,ref,subpxres) 
% finds which shift applied to 'box' gives to the best match to 'ref' at a subpixel resolution of 'subpxres'
% and returns a 4 element vector with the x y and z components of the shift and the value of the normalized cross-correlation. 
%
% box is expected to be smaller than ref and a shift of 0 corresponds to box and ref sharing the same center. 
% subpxres is the sub-pixel resolution required by the user (e.g. 1e-6) and the algorithm will proceed by dichotomies until reaching this precision level.
%
% Parts of this algorithm are inspired by: 
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).

if length(size(ref))~=3 || length(size(box))~=3
    error('The inputs are not 3D, terminating') % One possibility: if not 3D or if not a pixel of manoeuver available, add 1px of 0-padding on the ref, to avoid crashes. 
    return
end
correction=[0 0 0];
if any(floor(size(box)/2)==size(box)/2)
    correction=correction+double(floor(size(box)/2)==size(box)/2)/2;
    box=padarray(box,double(floor(size(box)/2)==size(box)/2),'post');
end
if any(floor(size(ref)/2)==size(ref)/2)
    correction=correction-double(floor(size(ref)/2)==size(ref)/2)/2;
    ref=padarray(ref,double(floor(size(ref)/2)==size(ref)/2),'post');
end
box_halfsize=floor(size(box)/2);
if ~exist('usfac','var')
    usfac = 1;
end

% The limit of displacement is implicit from the box sizes (we don't want
% any circular displacement, the reference is supposed to be wide enough to
% move the box within)
limit = floor(size(ref)/2)-floor(size(box)/2);

box=double(box);
ref=double(ref);

% Normalize the moving box:
N=numel(box);
box=(box-mean(box(:)));
tmp=mean(box(:).*box(:));
if tmp==0
    disp('Standard deviation of the moving box is 0. xcorr not defined (set to 0)');
    shiftxyz_CC=[0 0 0 0];
    return
else
    box=box/(sqrt(tmp));
end

% Compute the standard deviations of the reference under each box position:

% %%%%%% Naive algorithm:
% center_pos=ceil(size(ref)/2);
% ref_std=ones(size(ref));
% for i=-limit(1):limit(1)
%     for j=-limit(2):limit(2)
%         for k=-limit(3):limit(3)
%             tmp=ref(center_pos(1)-i-box_halfsize(1):center_pos(1)-i+box_halfsize(1),...
%                          center_pos(2)-j-box_halfsize(2):center_pos(2)-j+box_halfsize(2),...
%                          center_pos(3)-k-box_halfsize(3):center_pos(3)-k+box_halfsize(3));
%             ref_std(center_pos(1)+i,center_pos(2)+j,center_pos(3)+k)...
%                 =std(tmp(:));
%         end
%     end
% end

%%%%%%% Fast algorithm:
% Construct integral tables:
ref2=ref.*ref;
tmp=size(ref)+[1 1 1];
s=zeros(tmp);
s2=zeros(tmp);
for u=1:size(ref,1)
    for v=1:size(ref,2)
        for w=1:size(ref,3)
            s(u+1,v+1,w+1)=ref(u,v,w)+s(u,v+1,w+1)+s(u+1,v,w+1)+s(u+1,v+1,w)-s(u,v,w+1)-s(u,v+1,w)-s(u+1,v,w)+s(u,v,w);
            s2(u+1,v+1,w+1)=ref2(u,v,w)+s2(u,v+1,w+1)+s2(u+1,v,w+1)+s2(u+1,v+1,w)-s2(u,v,w+1)-s2(u,v+1,w)-s2(u+1,v,w)+s2(u,v,w);
        end
    end
end

ref_std=ones(size(ref));
center_pos=ceil(size(ref)/2);
for i=-limit(1):limit(1)
    for j=-limit(2):limit(2)
        for k=-limit(3):limit(3)
            xmin=center_pos(1)-i-box_halfsize(1);
            xmax=center_pos(1)-i+box_halfsize(1)+1;
            ymin=center_pos(2)-j-box_halfsize(2);
            ymax=center_pos(2)-j+box_halfsize(2)+1;
            zmin=center_pos(3)-k-box_halfsize(3);
            zmax=center_pos(3)-k+box_halfsize(3)+1;
            ref_sum=s(xmax,ymax,zmax)-s(xmin,ymax,zmax)-s(xmax,ymin,zmax)-s(xmax,ymax,zmin)+s(xmin,ymin,zmax)+s(xmax,ymin,zmin)+s(xmin,ymax,zmin)-s(xmin,ymin,zmin);
            squ_sum=s2(xmax,ymax,zmax)-s2(xmin,ymax,zmax)-s2(xmax,ymin,zmax)-s2(xmax,ymax,zmin)+s2(xmin,ymin,zmax)+s2(xmax,ymin,zmin)+s2(xmin,ymax,zmin)-s2(xmin,ymin,zmin);
            ref_mean=ref_sum/N;
            % If there is no contrast (std=0), CC is not defined: set CC to
            % 0 by setting the std to +Inf
            tmp=sqrt(squ_sum/N-ref_mean^2);
            if tmp==0
                ref_std(center_pos(1)+i,center_pos(2)+j,center_pos(3)+k)=+Inf;
            else
                ref_std(center_pos(1)+i,center_pos(2)+j,center_pos(3)+k)=tmp;
            end
        end
    end
end

% Make a circular permutation of the ref_std so that the std for no displacement is in 1,1,1:
ref_std(:,:,:)=ref_std([center_pos(1):end 1:center_pos(1)-1],...
                       [center_pos(2):end 1:center_pos(2)-1],...
                       [center_pos(3):end 1:center_pos(3)-1]);

% Pad the moving box for it to be the same size as the ref:
padsize=(size(ref)-size(box))/2;    %Expecting the size difference to be a multiple of 2 (same center)
if padsize~=floor(padsize)
    disp('The reference and box should have a size difference multiple of 2!');
    return
end
box=padarray(box,padsize);

% Compute the Fourrier transforms:
fft_box=fftn(box);
fft_box_ref=fftn(ref);

[nr,nc,nl]=size(fft_box_ref);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
Nl = ifftshift(-fix(nl/2):ceil(nl/2)-1);

% Single pixel accuracy registration
CC = ifftn(fft_box.*conj(fft_box_ref));
%CC=real(CC);    % Should be unnecessary as CC should be real, but just to be sure with the roundings.
CCnorm = CC./ref_std/N;     % Normalize with the standard deviation and the number of elements to get values between -1 and +1.
CCnorm(limit(1)+2:end-limit(1),:,:)=0;
CCnorm(:,limit(2)+2:end-limit(2),:)=0;
CCnorm(:,:,limit(3)+2:end-limit(3))=0;
[CCmax, ii] = max(CCnorm(:));
[row_shift_abs, col_shift_abs, lay_shift_abs] = ind2sub(size(CCnorm),ii);

% Now change shifts so that they represent relative shifts and not indices
row_shift_rel = Nr(row_shift_abs);
col_shift_rel = Nc(col_shift_abs);
lay_shift_rel = Nl(lay_shift_abs);

% Continue only if looking for subpixel resolution:
if subpxres<1
    %%% Now perform the subpixel accuracy calculation in the neighborhood using
    %%% the theoretical exact interpolation expression:
    
    % Make sure the subpixel doesn't get out of the search range:
    % (For a max displacement, we want the shift to be max-1 and the subpixel shift to be 1)
    if row_shift_rel==limit(1)
    row_shift_rel=row_shift_rel-1;
    row_shift_abs=row_shift_abs-1;
    end
    if col_shift_rel==limit(2)
        col_shift_rel=col_shift_rel-1;
        col_shift_abs=col_shift_abs-1;
    end
    if lay_shift_rel==limit(3)
        lay_shift_rel=lay_shift_rel-1;
        lay_shift_abs=lay_shift_abs-1;
    end
    
    row_shift_abs_save=row_shift_abs;
    col_shift_abs_save=col_shift_abs;
    lay_shift_abs_save=lay_shift_abs;
    row_shift_rel_save=row_shift_rel;
    col_shift_rel_save=col_shift_rel;
    lay_shift_rel_save=lay_shift_rel;
    
    % Look on both sides of the best integer value found and store the
    % cross-correlations found:
    CCsubpx=zeros([2,2,2]);
    xshifttmp=zeros([2,2,2]);
    yshifttmp=zeros([2,2,2]);
    zshifttmp=zeros([2,2,2]);
    
    for xreverse=[0,1]      % xreverse is 1 when looking for the subpixel optimum between -1 and 0 instead of 0 and 1. 
        for yreverse=[0,1]
            for zreverse=[0,1]
                
                if xreverse==1 && row_shift_rel_save>-limit(1)   % Only look on the backwards in each dimension if it doesn't get out of the search range.
                    row_shift_rel=row_shift_rel_save-1;
                    row_shift_abs=row_shift_abs_save-1;
                    if row_shift_abs==0
                        row_shift_abs=n(1);
                    end
                elseif xreverse==0
                    row_shift_rel=row_shift_rel_save;
                    row_shift_abs=row_shift_abs_save;
                end
                if yreverse==1 && col_shift_rel_save>-limit(2)
                    col_shift_rel=col_shift_rel_save-1;
                    col_shift_abs=col_shift_abs_save-1;
                    if col_shift_abs==0
                        col_shift_abs=n(2);
                    end
                elseif yreverse==0
                    col_shift_rel=col_shift_rel_save;
                    col_shift_abs=col_shift_abs_save;
                end
                if zreverse==1 && lay_shift_rel_save>-limit(3)
                    lay_shift_rel=lay_shift_rel_save-1;
                    lay_shift_abs=lay_shift_abs_save-1;
                    if lay_shift_abs==0
                        lay_shift_abs=n(3);
                    end
                elseif zreverse==0
                    lay_shift_rel=lay_shift_rel_save;
                    lay_shift_abs=lay_shift_abs_save;
                end
                
                % Pick up the already computed terms: 
                n=size(ref_std);
                if row_shift_abs==n(1)
                    row_indices=[n(1) 1];
                else
                    row_indices=[row_shift_abs row_shift_abs+1];
                end

                if col_shift_abs==n(2)
                    col_indices=[n(2) 1];
                else
                    col_indices=[col_shift_abs col_shift_abs+1];
                end

                if lay_shift_abs==n(3)
                    lay_indices=[n(3) 1];
                else
                    lay_indices=[lay_shift_abs lay_shift_abs+1];
                end

                % Known normalizations of the ref in the neighborhood:
                sigma_g=ref_std(row_indices,col_indices,lay_indices);

                % Known cross-correlations between the box and the ref in the neighborhood:
                cross_terms_fg=CC(row_indices,col_indices,lay_indices)/N;

                % Compute the cross-correlations between the ref and its neighbours in
                % the ref. Fill in also the squared standard deviations in the diagonal
                % (they are cross-correlations with self). 
                % Organization of the cross term matrix: 8x8 matrix with two indices, which
                % correspond to relative row, column, layer (i,j,k) of the two ref positions being
                % considered for the vector product:
                % index i j k
                % 1     0 0 0
                % 2     0 0 1
                % 3     0 1 0
                % 4     0 1 1
                % 5     1 0 0
                % 6     1 0 1
                % 7     1 1 0
                % 8     1 1 1
                % We have i=floor(mod(index-1,8)/4)=floor((index-1)/4);
                %         j=floor(mod(index-1,4)/2);
                %         k=floor(mod(index-1,2)/1)=mod(index-1,2);
                % And index=i*4+j*2+k+1;
                % The full matrix would be symmetrical, so only fill and use the upper triangle. 

                cross_terms_gg=zeros(8,8);

                box_ref_optimum=ref(center_pos(1)-row_shift_rel-box_halfsize(1)-1:center_pos(1)-row_shift_rel+box_halfsize(1),...
                                        center_pos(2)-col_shift_rel-box_halfsize(2)-1:center_pos(2)-col_shift_rel+box_halfsize(2),...
                                        center_pos(3)-lay_shift_rel-box_halfsize(3)-1:center_pos(3)-lay_shift_rel+box_halfsize(3));
                num=prod(size(box_ref_optimum)-1);

                for index1=1:8
                    i1=floor((index1-1)/4);
                    j1=floor(mod(index1-1,4)/2);
                    k1=mod(index1-1,2);
                    tmp1=box_ref_optimum(2-i1:end-i1,2-j1:end-j1,2-k1:end-k1);
                    tmp1=tmp1-mean(tmp1(:));
                    for index2=index1:8
                        if index1==index2
                            cross_terms_gg(index1,index2)=sigma_g(1+i1,1+j1,1+k1)^2;
                        else
                            i2=floor((index2-1)/4);
                            j2=floor(mod(index2-1,4)/2);
                            k2=mod(index2-1,2);
                            tmp2=box_ref_optimum(2-i2:end-i2,2-j2:end-j2,2-k2:end-k2);
                            tmp2=tmp2-mean(tmp2(:));
                            tmp3=tmp1.*tmp2;
                            cross_terms_gg(index1,index2)=sum(tmp3(:)/num);
                        end
                    end
                end

                % Iterative zoom around the best cross-correlation found
                % (~dichotomies method).
                step=1;
                x_best=0.5;
                y_best=0.5;
                z_best=0.5;
                while step>subpxres
                    step_old=step;
                    step=step/4;
                    
                    % Compute the subpixel symbolic expression:
                    x=x_best-step_old/2:step:x_best+step_old/2;
                    y=y_best-step_old/2:step:y_best+step_old/2;
                    z=z_best-step_old/2:step:z_best+step_old/2;
                    
                    [X,Y,Z]=meshgrid(x,y,z);
                    
                    L_sq=zeros(size(X));
                    for index1=1:8
                        i1=floor((index1-1)/4);
                        j1=floor(mod(index1-1,4)/2);
                        k1=mod(index1-1,2);
                        % Term with index2=index1:
                        L_sq=L_sq+cross_terms_gg(index1,index1)*(1-X).^(2-2*i1).*(1-Y).^(2-2*j1).*(1-Z).^(2-2*k1).*X.^(2*i1).*Y.^(2*j1).*Z.^(2*k1);
                        % Then cross-terms:
                        for index2=index1+1:8
                            i2=floor((index2-1)/4);
                            j2=floor(mod(index2-1,4)/2);
                            k2=mod(index2-1,2);
                            L_sq=L_sq+2*cross_terms_gg(index1,index2)*(1-X).^(2-i1-i2).*(1-Y).^(2-j1-j2).*(1-Z).^(2-k1-k2).*X.^(i1+i2).*Y.^(j1+j2).*Z.^(k1+k2);
                        end
                    end
                    L=sqrt(L_sq);
                    U=(1-X).*(1-Y).*(1-Z)*cross_terms_fg(1,1,1) ...
                            +   X .*(1-Y).*(1-Z)*cross_terms_fg(2,1,1) ...
                            +(1-X).*   Y .*(1-Z)*cross_terms_fg(1,2,1) ...
                            +   X .*   Y .*(1-Z)*cross_terms_fg(2,2,1) ...
                            +(1-X).*(1-Y).*   Z *cross_terms_fg(1,1,2) ...
                            +   X .*(1-Y).*   Z *cross_terms_fg(2,1,2) ...
                            +(1-X).*   Y .*   Z *cross_terms_fg(1,2,2) ...
                            +   X .*   Y .*   Z *cross_terms_fg(2,2,2);
                    CCxyz=U./L;
                    [CCmax, ii]=max(CCxyz(:));
                    x_best=X(ii);
                    y_best=Y(ii);
                    z_best=Z(ii);
                end
                
                CCsubpx(1+xreverse,1+yreverse,1+zreverse)=CCmax;
                xshifttmp(1+xreverse,1+yreverse,1+zreverse)=X(ii)-xreverse;
                yshifttmp(1+xreverse,1+yreverse,1+zreverse)=Y(ii)-yreverse;
                zshifttmp(1+xreverse,1+yreverse,1+zreverse)=Z(ii)-zreverse;
            end
        end
    end
    [CCmax,ii]=max(CCsubpx(:));
    xshift=xshifttmp(ii);
    yshift=yshifttmp(ii);
    zshift=zshifttmp(ii);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % If its only one row, column or layer, the shift along that dimension has no meaning. Set to zero.
    % if nr == 1,
    %     row_shift = 0;
    % end
    % if nc == 1,
    %     col_shift = 0;
    % end
    % if nl == 1,
    %     lay_shift = 0;
    % end

    shiftxyz_CC=[row_shift_rel_save+xshift,col_shift_rel_save+yshift,lay_shift_rel_save+zshift,CCmax]+[correction,0];
else
    shiftxyz_CC=[row_shift_rel,col_shift_rel,lay_shift_rel,CCmax]+[correction,0];
end
return


