disp('Opening the picture');
tic
if ~exist(input_stack) && exist(['data\' input_stack])
    input_stack=['data\' input_stack];
end
    
raw=bfopen(input_stack);    %If crashing here, be sure to transport the file to open to the current working folder, or change the address. 

% Put the images in a 3D matrix:
% File organization: XYCZT (5D hyperstack). 
experiment=1;
n=size(raw{experiment,1},1); % n is the number of 2D planes (n_layers x n_timepoints x n_colors).
n_px=size(raw{experiment,1}{1,1});
n_layers=n/n_timepoints/n_colors;

% Probably not needed anymore:
% clear img
% img=zeros(n_px(1),n_px(2),n,class(raw{experiment,1}{1,1}));
% for i=1:n
%     img(:,:,i)=raw{experiment,1}{i,1};
% end
% img=double(img(:,:,tracking_channel:n_colors:end));

% img2 has the original tiff organization (e.g. red fibers and
% green cell), used for the display functions:
img=zeros(n_px(1),n_px(2),n_colors,n_layers,n_timepoints,'uint8');
for t=1:n_timepoints
    for z=1:n_layers
        for c=1:n_colors
            ii=n_layers*n_colors*(t-1)+n_colors*(z-1)+c;
            img(:,:,c,z,t)=raw{experiment,1}{ii,1}; % img has indices: x, y, c(channel), z(layer), t(timepoint)
        end
    end
end

clear raw

toc
