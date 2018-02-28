function [ output_args ] = simple_quiver(U,varargin)
% simple_quiver(U) or (U,'color'), does a 3D quiver plot of a vector field
% U(Px,Py,Pz,Vcomponent) with Px, Py, Pz indices of the arrow position and 
% Vcomponent a value between 1 and 3 pointing to the components (xyz) of the
% vectors. 

s=size(U);
if length(s)>4
    U=squeeze(U);
end

s=size(U);
[yy,xx,zz]=meshgrid(1:s(1),1:s(2),1:s(3));
if length(varargin)==1
    quiver3(xx,yy,zz,U(:,:,:,1),U(:,:,:,2),U(:,:,:,3),varargin{1});
else
    quiver3(xx,yy,zz,U(:,:,:,1),U(:,:,:,2),U(:,:,:,3));
end

end

