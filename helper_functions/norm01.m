function [M_normalized] = norm01(M)
%Normalizes the matrix so its values lie between 0 and 1. 

M=double(M);
M_normalized=(M-min(M(:)))/(max(M(:))-min(M(:)));

end
