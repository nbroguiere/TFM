function [ output_args ] = imsave3D( image3D, name )
%Export a 3D stack bw or color to the harddrive. 
%3D stack is in the form of a 3-entry matrix. Values between 0 and 1. 
%Nicolas Broguiere 22.04.2014

if(exist(name,'file'))
    delete(name);
end
if length(size(image3D))==3  %bw xyz or xyt
    for i=1:size(image3D,3)
        imwrite(image3D(:,:,i), name, 'WriteMode','append');
    end
elseif length(size(image3D))==4  %color xycz or xyct
    if size(image3D,3)==3   % Normal RGB picture
        for i=1:size(image3D,4)
            imwrite(image3D(:,:,:,i), name, 'WriteMode','append');
        end
    else                    % Composite image not RGB: separate the channels in alternate layers.
        image3D=reshape(image3D,size(image3D,1),size(image3D,2),size(image3D,3)*size(image3D,4));
        for i=1:size(image3D,3)
            imwrite(image3D(:,:,i), name, 'WriteMode','append');
        end
    end
elseif length(size(image3D))==5  % color timelapse xyczt
    if size(image3D,3)==3   % Normal RGB picture
        for t=1:size(image3D,5)
            for z=1:size(image3D,4)
                imwrite(image3D(:,:,:,z,t), name, 'WriteMode','append');
            end
        end
    else                    % Composite image not RGB: separate the channels in alternate layers.
        image3D=reshape(image3D,size(image3D,1),size(image3D,2),size(image3D,3)*size(image3D,4)*size(image3D,5));
        for i=1:size(image3D,3)
            imwrite(image3D(:,:,i), name, 'WriteMode','append');
        end
    end
end
end

function [ output_args ] = imsave3D( image3D, name )
%Export a 3D stack bw or color to the harddrive. 
%3D stack is in the form of a 3-entry matrix. Values between 0 and 1. 
%Nicolas Broguiere 22.04.2014

if(exist(name,'file'))
    delete(name);
end
if length(size(image3D))==3  %bw xyz or xyt
    for i=1:size(image3D,3)
        imwrite(image3D(:,:,i), name, 'WriteMode','append');
    end
elseif length(size(image3D))==4  %color xycz or xyct
    if size(image3D,3)==3   % Normal RGB picture
        for i=1:size(image3D,4)
            imwrite(image3D(:,:,:,i), name, 'WriteMode','append');
        end
    else                    % Composite image not RGB: separate the channels in alternate layers.
        image3D=reshape(image3D,size(image3D,1),size(image3D,2),size(image3D,3)*size(image3D,4));
        for i=1:size(image3D,3)
            imwrite(image3D(:,:,i), name, 'WriteMode','append');
        end
    end
elseif length(size(image3D))==5  % color timelapse xyczt
    if size(image3D,3)==3   % Normal RGB picture
        for t=1:size(image3D,5)
            for z=1:size(image3D,4)
                imwrite(image3D(:,:,:,z,t), name, 'WriteMode','append');
            end
        end
    else                    % Composite image not RGB: separate the channels in alternate layers.
        image3D=reshape(image3D,size(image3D,1),size(image3D,2),size(image3D,3)*size(image3D,4)*size(image3D,5));
        for i=1:size(image3D,3)
            imwrite(image3D(:,:,i), name, 'WriteMode','append');
        end
    end
end
end

