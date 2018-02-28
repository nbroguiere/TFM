function [output] = xcorr3subpxnormQuadInt(box,ref)
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
    output=[0 0 0 0];
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

%%% Now perform the subpixel accuracy calculation in the neighborhood using
%%% interpolation over neighbouring pixels with a second order polynomial:
% The max is xmax=(v-1-v1)/2(v1+v-1-2v0)

% Make sure to handle the limit cases properly: dont look for the subpixel
% shift if the integer pixel shift is already at the limit, and when
% reaching the end of the matrix go back to first position and vice versa. 

n=size(ref_std);
if abs(row_shift_rel)==limit(1)
    xshift=0;
else
    if row_shift_abs==n(1)
        row_indices=[n(1)-1 n(1) 1];
    elseif row_shift_abs==1
        row_indices=[n(1) 1 2];
    else
        row_indices=[row_shift_abs-1 row_shift_abs row_shift_abs+1];
    end
    CCx=CCnorm(row_indices(1:3),col_shift_abs,lay_shift_abs); % CCx contains cross-correlation in previous, middle, and next pixels in this order
    xshift=(CCx(1)-CCx(3))/(2*(CCx(3)+CCx(1)-2*CCx(2)));
end
if abs(col_shift_rel)==limit(2)
    yshift=0;
else
    if col_shift_abs==n(2)
        col_indices=[n(2)-1 n(2) 1];
    elseif col_shift_abs==1
        col_indices=[n(2) 1 2];
    else
        col_indices=[col_shift_abs-1 col_shift_abs col_shift_abs+1];
    end
    CCy=CCnorm(row_shift_abs,col_indices(1:3),lay_shift_abs); % CCy contains cross-correlation in previous, middle, and next pixels in this order
    yshift=(CCy(1)-CCy(3))/(2*(CCy(3)+CCy(1)-2*CCy(2)));
end
if abs(lay_shift_rel)==limit(3)
    zshift=0;
else
    if lay_shift_abs==n(3)
        lay_indices=[n(3)-1 n(3) 1];
    elseif lay_shift_abs==1
        lay_indices=[n(3) 1 2];
    else
        lay_indices=[lay_shift_abs-1 lay_shift_abs lay_shift_abs+1];
    end
    CCz=CCnorm(row_shift_abs,col_shift_abs,lay_indices(1:3)); % CCz contains cross-correlation in previous, middle, and next pixels in this order
    zshift=(CCz(1)-CCz(3))/(2*(CCz(3)+CCz(1)-2*CCz(2)));
end

output=[row_shift_rel+xshift,col_shift_rel+yshift,lay_shift_rel+zshift,CCmax]+[correction,0];
return


