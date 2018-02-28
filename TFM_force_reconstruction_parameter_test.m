tic
% Choice of the testing model:
model=3;             % Model 1: Half plane under compression. Model 2: Half plane under shear. Model 3: Sphere under compression. 
direction=1;         % Compression or shear in the x, y or z direction (1,2,3 respectively). 
strain_applied=10;   % in percentage (of the diameter of the sphere or thickness of the half plane)
r0=15;               % for the sphere model only, half-diameter of the sphere, in um
noise_level=100;      % in percentage of the strain. (for example 20% noise: adding from -10% of max strain to +10% of max strain with a constant probability distribution). 

% Displacement field fitting parameters:
alpha_F_cell_surface=1e-9;      % Relative importance of minimizing the forces on the cell surface 
alpha_F_matrix=1e9;             % Relative importance of minimizing the forces compared to fitting the data in the matrix. Important note: in trials, when going up to 1e14, all fine. When going to 1e15, very different solution. Probably problems with the roundings. 
alpha_smoothing=100;            % Relative importance of minimizing on the cell surface the difference between the force on a position and the average force on neighboring positions.
alpha_z=1;                      % Relative importance of good fitting on z values (which can be less accurate than xy values). Values between 0 (ignore z values) and 1 (put the same penalty as xy values if the fitting is not good).

% Mechanical parameters (Lamé coefficients):
lambda=0;   % in N/m^2 = pN/um^2
mu=100;     % in N/m^2 = pN/um^2

step_size_xy=4;                % With current parameters, starts to gets slow under 7 px
step_size_z=step_size_xy;

pxsize_xy=0.2;        % um/px
pxsize_z=pxsize_xy;   % um
area_size=40;         % um, size of the complete working field. 
n_px=[floor(area_size/pxsize_xy) floor(area_size/pxsize_xy)];
n_layers=floor(area_size/pxsize_z);
sigma=0.4;            % um, standard deviation of the 3D Gaussian smoothing applied before meshing the sphere surface (in model 3 only). 
patch_reduction=0.01; % Fraction of the facets kept (before reduction, the facet numberis close to the number of pixels on the surface). 

dx=step_size_xy*pxsize_xy;
dy=step_size_xy*pxsize_xy;
dz=step_size_z*pxsize_z;

% Adimensionalize the coefficients so they don't need to be changed for different step size etc. 
F0=mu/(dx*dy*dz)^(1/3);                                 % Typical volumic force used for adimensionalization.
alpha_F_cell_surface=alpha_F_cell_surface/F0^2;         % Relative importance of minimizing the forces on the cell surface 
alpha_F_matrix=alpha_F_matrix/F0^2;                     % Relative importance of minimizing the forces compared to fitting the data in the matrix. Important note: in trials, when going up to 1e14, all fine. When going to 1e15, very different solution. Probably problems with the roundings. 
alpha_smoothing=alpha_smoothing/F0^2;                   % Relative importance of minimizing on the cell surface the difference between the force on a position and the average force on neighboring positions.

% Location of the scripts:
addpath(genpath('C:\Users\nicki\Documents\PhD DZNE\Traction force microscopy\MATLAB functions and plugins'));

xrange=1:step_size_xy:n_px(1);
yrange=1:step_size_xy:n_px(2);
zrange=1:step_size_z:n_layers;
nx=length(xrange);
ny=length(yrange);
nz=length(zrange);

[Py, Px, Pz]=meshgrid(xrange,yrange,zrange);
[yy, xx, zz]=meshgrid(1:n_px(1),1:n_px(2),1:n_layers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define the simulated displacement fields for the various models %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% define pushing on half of the plane (20 um length):
if model==1
    U0=area_size*strain_applied/100;% in um. 
    if direction==1
        img_threshold=xx>=n_px(1)/2; 
        Ux=U0*Px/n_px(1);   % in um. 
        Uy=zeros(size(Px)); % in um. 
        Uz=zeros(size(Px)); % in um. 
    elseif direction==2
        img_threshold=yy>=n_px(2)/2; 
        Ux=zeros(size(Px));   % in um. 
        Uy=U0*Py/n_px(2); % in um. 
        Uz=zeros(size(Px)); % in um. 
    elseif direction==3
        img_threshold=zz>=n_layers/2; 
        Ux=zeros(size(Px));   % in um. 
        Uy=zeros(size(Px)); % in um. 
        Uz=U0*Pz/n_layers; % in um. 
    end
    
    % Area of the cell surface:
    S=area_size*area_size;
end
%%%%%%%% define shearing half of the plane (20 um length):
if model==2
    U0=area_size*strain_applied/100;    % in um. 
    if direction==1
        img_threshold=xx>=n_px(1)/2; 
        Ux=zeros(size(Px)); % in um. 
        Uy=U0*Px/n_px(1);   % in um. 
        Uz=zeros(size(Px)); % in um. 
    elseif direction==2
        img_threshold=yy>=n_px(2)/2; 
        Ux=zeros(size(Px));   % in um. 
        Uy=zeros(size(Px)); % in um. 
        Uz=U0*Py/n_px(2); % in um. 
    elseif direction==3
        img_threshold=zz>=n_layers/2; 
        Ux=U0*Pz/n_layers; % in um. 
        Uy=zeros(size(Px)); % in um. 
        Uz=zeros(size(Px));   % in um. 
    end
    % Area of the cell surface:
    S=area_size*area_size;
end
%%%%%%%% define spherical compression:
if model==3
    U0=r0*strain_applied/100;   % in um
    
    img_threshold=(pxsize_xy^2*(xx-n_px(1)/2).^2+pxsize_xy^2*(yy-n_px(2)/2).^2+pxsize_z^2*(zz-n_layers/2).^2)<r0.^2;
    
    % T0*r0^2/r^2*e_r:
    r2=pxsize_xy^2*(Px-n_px(1)/2).^2+pxsize_xy^2*(Py-n_px(2)/2).^2+pxsize_z^2*(Pz-n_layers/2).^2;
    r3=sqrt(r2).*r2;
    Ux=U0*r0^2*pxsize_xy*(Px-n_px(1)/2)./r3; % in um. 
    Uy=U0*r0^2*pxsize_xy*(Py-n_px(2)/2)./r3; % in um. 
    Uz=U0*r0^2*pxsize_z*(Pz-n_layers/2)./r3; % in um. 
    
    % Area of the cell surface:
    S=4*pi*r0^2;
end

% Save the displacements without noise (=theoretical solution):
Ux_theory=Ux;
Uy_theory=Uy;
Uz_theory=Uz;

% Add noise:
noise_um=noise_level/100*strain_applied/100*r0;
Ux=Ux+(rand(size(Ux))-0.5)*noise_um;
Uy=Uy+(rand(size(Uy))-0.5)*noise_um;
Uz=Uz+(rand(size(Uz))-0.5)*noise_um;

%%%%%%%% end of the definition

%%%%%% The cell is defined as all the positions containing some positive
%%%%%% values in the img_threshold. 
%%%%%% The surface of the cells is found from the cell position as the
%%%%%% cell positions which have a neighbor outside the cell. 
%%%%%% The cell position is finally refined to exclude the cell surface. 

%disp('finding cell surface')
cell_position=zeros(size(Px),'logical');    % indices: i j k t
cell_surface=zeros(size(Px),'logical');    % indices: i j k t
sxyp=floor(step_size_xy/2);
sxym=floor((step_size_xy-1)/2);
szp=floor(step_size_z/2);
szm=floor((step_size_z-1)/2);
numelbox=step_size_xy.^2*step_size_z;
neighborhood=zeros([3,3,3],'logical');
neighborhood(1:3,2,2)=1;
neighborhood(2,1:3,2)=1;
neighborhood(2,2,1:3)=1;
for k=1:nz
    z=zrange(k);
    zmin=max(z-szm,1);
    zmax=min(z+szp,n_layers);
    for j=1:ny
        y=yrange(j);
        ymin=max(y-sxym,1);
        ymax=min(y+sxyp,n_px(2));
        for i=1:nx
            x=xrange(i);
            xmin=max(x-sxym,1);
            xmax=min(x+sxyp,n_px(1));

            tmp=img_threshold(xmin:xmax,ymin:ymax,zmin:zmax);
            tmp2=sum(double(tmp(:)));
            if tmp2>0 %&& tmp2<numelbox
                cell_position(i,j,k)=1;
            end
        end
    end
end
for i=1:nx
    for j=1:ny
        for k=1:nz
            if cell_position(i,j,k)
                tmp=cell_position(max(i-1,1):min(i+1,nx),max(j-1,1):min(j+1,ny),max(k-1,1):min(k+1,nz));
                i1=1;i2=3;j1=1;j2=3;k1=1;k2=3;
                if i==1
                    i1=2;
                end
                if i==nx
                    i2=2;
                end
                if j==1
                    j1=2;
                end
                if j==ny
                    j2=2;
                end
                if k==1
                    k1=2;
                end
                if k==nz
                    k2=2;
                end
                tmp2=neighborhood(i1:i2,j1:j2,k1:k2);
                tmp3=tmp.*tmp2;
                if sum(tmp3(:))<sum(tmp2(:))
                    cell_surface(i,j,k)=true;
                end
            end
        end
    end
end
cell_position=xor(cell_position,cell_surface);

Ux(cell_position)=nan;
Uy(cell_position)=nan;
Uz(cell_position)=nan;

Ux_theory(cell_position)=nan;
Uy_theory(cell_position)=nan;
Uz_theory(cell_position)=nan;
U_theory=cat(4,Ux_theory,Uy_theory,Uz_theory);

%disp('Building the matrix M of the mechanical model so that Force=M*Displacement');
%%%%%%%%% Build the matrix giving the force from the displacement (in vector form):
% Generate template matrices to build the blocks: 
% (put the neighboring value in the border conditions: equivalent to no deformations outside of the field of analysis):
% (one way to see it is outside of the simulations, only rigid bodies.
% Another way to see it is that mathematically with a solid border, the
% second derivative 1 -2 1 comes down to a first derivative, -1 1, so we
% find the stress rather than volumic forces). 
n_block=nx*ny*nz;
Identity=spdiags(ones(n_block,1),0,n_block,n_block);
next=1;
prev=-1;
e=ones(n_block-1,1);
e(nx:nx:n_block)=0;
Mxn=spdiags([0; e],next,n_block,n_block);   % Mxn is the matrix that picks up the next element in the x direction
Mxp=spdiags([e; 0],prev,n_block,n_block);   % Mxp is the matrix that picks up the prev element in the x direction
for j=0:ny-1
    for k=0:nz-1
        index=1+j*nx+k*nx*ny;
        Mxp(index,index)=1;
        index=nx+j*nx+k*nx*ny;
        Mxn(index,index)=1;
    end
end
e1=ones(nx*(ny-1),1);
e2=zeros(nx,1);
e=repmat([e1;e2],nz,1);
e=e(1:end-nx);
Myn=spdiags([zeros(nx,1);e],next*nx,n_block,n_block);
Myp=spdiags([e;zeros(nx,1)],prev*nx,n_block,n_block);
for k=0:nz-1
    index=1+k*nx*ny;
    for i=1:nx
        Myp(index,index)=1;
        index=index+1;
    end
    index=1+nx*(ny-1)+k*nx*ny;
    for i=1:nx
        Myn(index,index)=1;
        index=index+1;
    end
end
e=ones(n_block,1);
Mzn=spdiags(e,next*nx*ny,n_block,n_block);
Mzp=spdiags(e,prev*nx*ny,n_block,n_block);
for ij=1:nx*ny
    Mzp(ij,ij)=1;
end
for ij=1+nx*ny*(nz-1):nx*ny*nz
    Mzn(ij,ij)=1;
end

% Compute the first derivative matrices (Before altering the matrices for reflected values on borders in second derivatives)
I=Identity;
% derivative on the right: 
Drx=(Mxn-I)/dx;
Dry=(Myn-I)/dy;
Drz=(Mzn-I)/dz;
% derivative on the left:
Dlx=(I-Mxp)/dx;
Dly=(I-Myp)/dy;
Dlz=(I-Mzp)/dz;

% Fix the border conditions on the surface of the cells to also have reflected values: 
indices=find(cell_surface);

%cell_incl_surface=or(cell_position,cell_surface);
for index=indices'
    [i,j,k]=ind2sub([nx,ny,nz],index);
    % Find the direction of the surface in x and fix what needs to be:
    if i==nx
    elseif cell_position(i+1,j,k)   % If the cell is just next, calls for the next value should return the local value. 
        Mxn(i+nx*(j-1)+nx*ny*(k-1),i+1+nx*(j-1)+nx*ny*(k-1))=0;
        Mxn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
    if i==1
    elseif cell_position(i-1,j,k)   % If the cell is just before, calls for the previous value should return the local value. 
        Mxp(i+nx*(j-1)+nx*ny*(k-1),i-1+nx*(j-1)+nx*ny*(k-1))=0;
        Mxp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
    % Same in y:
    if j==ny
    elseif cell_position(i,j+1,k)   % If the cell is just next, calls for the next value should return the local value. 
        Myn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j)+nx*ny*(k-1))=0;
        Myn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
    if j==1
    elseif cell_position(i,j-1,k)   % If the cell is just before, calls for the previous value should return the local value. 
        Myp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-2)+nx*ny*(k-1))=0;
        Myp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
    % Same in z:
    if k==nz
    elseif cell_position(i,j,k+1)   % If the cell is just next, calls for the next value should return the local value. 
        Mzn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k))=0;
        Mzn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
    if k==1
    elseif cell_position(i,j,k-1)   % If the cell is just before, calls for the previous value should return the local value. 
        Mzp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-2))=0;
        Mzp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
    end
end

% Then build the convolution filters h corresponding to the various second derivatives:
second_derivative_straight=[1 -2 1];
second_derivative_crossed=1/2*[1 -1 0; -1 2 -1; 0 -1 1];
h=cell(3,3);
h{1,1}=zeros(3,3,3);
h{1,1}(:,2,2)=second_derivative_straight/(dx*dx);
h{1,2}=zeros(3,3,3);
h{1,2}(:,:,2)=second_derivative_crossed/(dx*dy);
h{1,3}=zeros(3,3,3);
h{1,3}(:,2,:)=second_derivative_crossed/(dx*dz);
h{2,2}=zeros(3,3,3);
h{2,2}(2,:,2)=second_derivative_straight/(dy*dy);
h{2,3}=zeros(3,3,3);
h{2,3}(2,:,:)=second_derivative_crossed/(dy*dz);
h{3,3}=zeros(3,3,3);
h{3,3}(2,2,:)=second_derivative_straight/(dz*dz);

% Use the convolution filters and the displacement matrices to create the
% block matrices H corresponding to the various derivatives:
H=cell(3,3);
for i=1:3
    for j=i:3
        if i==j
            H{i,j}=h{i,j}(2,2,2)*Identity...
                  +h{i,j}(1,2,2)*Mxp+h{i,j}(3,2,2)*Mxn...
                  +h{i,j}(2,1,2)*Myp+h{i,j}(2,3,2)*Myn...
                  +h{i,j}(2,2,1)*Mzp+h{i,j}(2,2,3)*Mzn;
        else
            H{i,j}=h{i,j}(2,2,2)*Identity...
                  +h{i,j}(1,2,2)*Mxp+h{i,j}(3,2,2)*Mxn...
                  +h{i,j}(2,1,2)*Myp+h{i,j}(2,3,2)*Myn...
                  +h{i,j}(2,2,1)*Mzp+h{i,j}(2,2,3)*Mzn...
                  +h{i,j}(1,2,1)*Mxp*Mzp+h{i,j}(1,2,3)*Mxp*Mzn...
                  +h{i,j}(3,2,1)*Mxn*Mzp+h{i,j}(3,2,3)*Mxn*Mzn...
                  +h{i,j}(1,1,2)*Mxp*Myp+h{i,j}(1,3,2)*Mxp*Myn...
                  +h{i,j}(3,1,2)*Mxn*Myp+h{i,j}(3,3,2)*Mxn*Myn...
                  +h{i,j}(2,1,1)*Myp*Mzp+h{i,j}(2,1,3)*Myp*Mzn...
                  +h{i,j}(2,3,1)*Myn*Mzp+h{i,j}(2,3,3)*Myn*Mzn;
        end
    end
end
for i=2:3
    for j=1:i-1
        H{i,j}=H{j,i};
    end
end

% %%%%%% test
% tmp=H{2,2};
% H{2,2}=H{1,1};
% H{1,1}=tmp;
% %%%%%%%%

% Put them all together to build the nine blocks to assemble to get the final matrix M satisfying F=M*U:
% In this notation, i represents x y or z. 
% block{i,j} gives the component of Fi coming from Uj.
% Overall, the physical equation is 
%     eps=1/2*(du_i/dx_j+du_j/dx_i)
%     sigma=2*mu*eps+lambda*tr(eps)*I
%     F=-sum(dsigma_ij/dx_i)
block=cell(3,3);
for i=1:3
    for j=1:3
        if i==j
            block{i,j}=-(lambda+2*mu)*H{i,i};
            for k=1:3
                if k~=i
                    block{i,j}=block{i,j}-mu*H{k,k}; % lambda and mu are in pN/um^2 and H in um^-2 so M is in pN/um^4
                end
            end
        else
            block{i,j}=-(lambda+mu)*H{i,j};
        end
    end
end
M=[block{1,1} block{1,2} block{1,3}; block{2,1} block{2,2} block{2,3}; block{3,1} block{3,2} block{3,3}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/ building force matrix

%%%%%%%% Similarly, build the average-of-neighboring-forces matrix %%%%%%%%
% Define the standard convolution filter defining neighbours (here: 27),
% not normalized (make it later on in the circle to take into account the
% actual number of neighbours for each surface position line). 
%disp('Building the smoothing matrix');
h_avg=ones(3,3,3);
h_avg(2,2)=0;
H_avg=h_avg(2,2,2)*Identity...
     +h_avg(1,2,2)*Mxp+h_avg(3,2,2)*Mxn...
     +h_avg(2,1,2)*Myp+h_avg(2,3,2)*Myn...
     +h_avg(2,2,1)*Mzp+h_avg(2,2,3)*Mzn...
     +h_avg(1,2,1)*Mxp*Mzp+h_avg(1,2,3)*Mxp*Mzn...
     +h_avg(3,2,1)*Mxn*Mzp+h_avg(3,2,3)*Mxn*Mzn...
     +h_avg(1,1,2)*Mxp*Myp+h_avg(1,3,2)*Mxp*Myn...
     +h_avg(3,1,2)*Mxn*Myp+h_avg(3,3,2)*Mxn*Myn...
     +h_avg(2,1,1)*Myp*Mzp+h_avg(2,1,3)*Myp*Mzn...
     +h_avg(2,3,1)*Myn*Mzp+h_avg(2,3,3)*Myn*Mzn...
     +h_avg(1,1,1)*Mxp*Myp*Mzp+h_avg(1,3,1)*Mxp*Myn*Mzp...
     +h_avg(3,1,1)*Mxn*Myp*Mzp+h_avg(3,3,1)*Mxn*Myn*Mzp...
     +h_avg(1,1,3)*Mxp*Myp*Mzn+h_avg(1,3,3)*Mxp*Myn*Mzn...
     +h_avg(3,1,3)*Mxn*Myp*Mzn+h_avg(3,3,3)*Mxn*Myn*Mzn;

%disp('Finding the displacement field generated from forces restricted to the cell surface that matches the experimental displacement field best'); 
Ufit=zeros(nx,ny,nz,3);
Ffit=zeros(nx,ny,nz,3);
%for t=timepoints_to_process
    % Create U which is the field of displacement, but in 1D-storage form:
    Uexp=[Ux(:); Uy(:); Uz(:)]; % in um
    n_U=nx*ny*nz*3;
    
    % Create a weight matrix with zeroes when the cross-correlation
    % coefficient has been deemed unacceptable:
    acceptable=~isnan(Uexp);
    % (also apply a coefficient for how important is good fit of z data compared
    % to x and y data (note that z components of vectors are not as good
    % because of the quality of the raw data is far lesser):
    tmp=double(acceptable);
    tmp(n_block*2+1:n_block*3)=tmp(n_block*2+1:n_block*3)*sqrt(alpha_z);
    A=spdiags(tmp,0,sparse(n_U,n_U));
    
    % Create the coefficient matrix B for the force minimization (min of F'*B'*B*F=U'*M'*B'*B*M*U=norm(f(U)))
    % Note that B includes different coefficients for the matrix (alpha high, forces forbidden) and on the cell surface (alpha low, forces authorized)
    % For cell surface and borders authorize forces (excluding borders which are inside the cell):
    tmp=zeros(nx-2,ny-2,nz-2);
    borders=logical(padarray(tmp,[1,1,1],1));
    surface_and_borders=or(cell_surface,and(borders,~cell_position));   % cell surface and borders which are not in the cell.
    if alpha_F_cell_surface>=0
        tmp=surface_and_borders(:,:,:)*sqrt(alpha_F_cell_surface);
    else
        tmp=surface_and_borders(:,:,:)*sqrt(-alpha_F_cell_surface);
    end    
    tmp2=[tmp(:);tmp(:);tmp(:)];
    B1=spdiags(tmp2(:),0,sparse(n_U,n_U));
    % For matrix only:
    tmp=~surface_and_borders(:,:,:)*sqrt(alpha_F_matrix);
    tmp2=[tmp(:);tmp(:);tmp(:)];
    B2=spdiags(tmp2(:),0,sparse(n_U,n_U));
    B=B1+B2;
    
    % Create the smoothing matrix D for the force regularization on the cell surface (min of F'*C'*C*F=U'*M'*C'*C*M*U=norm(F-Favg))
    tmp=cell_surface(:,:,:);  % Pick up the cell surface positions at this timepoint.
    tmp=sparse(tmp(:)');                % Linearize the image to a vector
    H_avg_norm=tmp'.*H_avg.*tmp;           % Take a contribution coefficient of 0 for the positions not on the cell surface, and dont send values to positions outside the cell surface.
    H_avg_norm=H_avg_norm./max(sum(H_avg_norm,2),eps);  % Normalize each line so the average is actually an average (e.g. if there are 3 neighours, divide by 3). When there is a division by zero, divide by epsilon instead. 
    tmp=Identity-H_avg_norm;
    D=sqrt(alpha_smoothing)*blkdiag(tmp,tmp,tmp);

    %%%%% Deprecated comment and code, now border conditions are fixed. 
%     % On edges, the forces are not well defined (would need a value outside of the window
%     % to compute the derivatives). So make sure the equation is 0=0 so that the border values
%     % of Ufit are determined by the equations in the second row from the
%     % border, where forces are correctly defined. 
%     if ~force_null_border
%         tmp=ones(nx-2,ny-2,nz-2);
%         tmp=padarray(tmp,[1,1,1],eps);
%         tmp2=[tmp(:);tmp(:);tmp(:)];
%         no_constrain_on_border=spdiags(tmp2(:),0,sparse(n_U,n_U));
%         tmp=ones(nx-4,ny-4,nz-4);
%         tmp=padarray(tmp,[2,2,2],0);
%         tmp2=[tmp(:);tmp(:);tmp(:)];
%         no_constrain_on_border2=spdiags(tmp2(:),0,sparse(n_U,n_U));
%         A=A*no_constrain_on_border;     % importance of fitting the displacement value
%         B1=B1*no_constrain_on_border;   % importance of minimizing the forces in the matrix
%         B2=B2*no_constrain_on_border;   % importance of minimizing the forces on the cell surface
%         D=D*no_constrain_on_border;    % importance of smoothing the forces
%     end
    
    % Then remove the nans from the experimental data, the data wont be used but nans would propagate:
    Uexp(isnan(Uexp))=0;

    % Construct the matrices in a way that makes sure matlab recognizes
    % them as Hermitian for the linear solving: 
    tmp=A*Uexp;
    tmp2_1=B1*M;
    tmp2_2=B2*M;
    tmp3=D*M;
    if alpha_F_cell_surface>=0
        tmp4=A+tmp2_1'*tmp2_1+tmp2_2'*tmp2_2+tmp3'*tmp3;
    else
        tmp4=A-tmp2_1'*tmp2_1+tmp2_2'*tmp2_2+tmp3'*tmp3;      % This way of constructing the matrix makes sure matlab recognizes it as Hermitian, makes solving the system faster. 
    end
    Vfit=tmp4\tmp;      % Solving the finite element model
    % Faster but imprecise solver:
    %Vfit=cgs(tmp4,tmp,1e-15,30,[],[],tmp);
    Ffit_vector=M*Vfit;    % Computing the associated force field in pN/um^3
%     if ~force_null_border
%         Ufit(:,:,:,:)=reshape(no_constrain_on_border*Vfit,nx,ny,nz,3);  % Puting back from vector form to normal 3D vector field organization
%         Ffit(:,:,:,:)=reshape(no_constrain_on_border*Ffit_vector,nx,ny,nz,3);
%     else
        Ufit(:,:,:,:)=reshape(Vfit,nx,ny,nz,3);  % Puting back from vector form to normal 3D vector field organization
        Ffit(:,:,:,:)=reshape(Ffit_vector,nx,ny,nz,3);
%     end
%end
%clear tmpx tmpy tmpz Ffit_vector tmp tmp2 tmp3 tmp4 A sqrtB border_null
Uexp=reshape(Uexp,size(Ufit));  % in um
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Finding the surface of the cell  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth the picture (to have a better patch with the sphere):
if model==3
    %disp('Smoothing the segmented cell surface')
    img_threshold_smooth=double(img_threshold);
    box_range=-floor(3*sigma/pxsize_xy):floor(3*sigma/pxsize_xy);
    filter=normpdf(box_range,0,sigma/pxsize_xy);
    img_threshold_smooth=imfilter(img_threshold_smooth,filter);
    img_threshold_smooth=imfilter(img_threshold_smooth,filter');
    box_range=-floor(3*sigma/pxsize_z):floor(3*sigma/pxsize_z);
    filter=zeros(1,1,length(box_range));
    filter(:)=normpdf(box_range,0,sigma/pxsize_z);
    img_threshold_smooth=imfilter(img_threshold_smooth,filter);
end

[xx, yy, zz]=meshgrid((1:n_px(2))*pxsize_xy,(1:n_px(1))*pxsize_xy,(1:n_layers)*pxsize_z);
fv=isosurface(xx,yy,zz,img_threshold_smooth,0.5);
if model==3
    fv=reducepatch(fv,patch_reduction);
end
% Find the surface area of the cell (i.e. the patch) and the facet centers and normal unit vectors:
a = fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,1),:);
b = fv.vertices(fv.faces(:,3),:)-fv.vertices(fv.faces(:,1),:);
c = cross(a, b, 2);
area = 1/2 * sum(sqrt(sum(c.^2, 2)));
%disp(['The surface area is ' num2str(area) ' um^2.'])
tmp=sqrt(sum(c.^2,2));
normals=-c./cat(2,tmp,tmp,tmp);
facet_centers=(fv.vertices(fv.faces(:,1),:)+fv.vertices(fv.faces(:,2),:)+fv.vertices(fv.faces(:,3),:))/3;   % coordinates in um. 
% to check the normals, use:
%quiver3(facet_centers(:,1),facet_centers(:,2),facet_centers(:,3),normals(:,1),normals(:,2),normals(:,3),2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=cat(4,cell_position,cell_position,cell_position);
Ufit(tmp)=nan;
Uexp(tmp)=nan;

% F_avg=blkdiag(H_avg_norm,H_avg_norm,H_avg_norm)*Ffit(:);
% F_avg=reshape(F_avg,size(Ffit));

Ufit_norm=sqrt(Ufit(:,:,:,1).^2+Ufit(:,:,:,2).^2+Ufit(:,:,:,3).^2); % Ufit is in same units as T, so um here. 
Uexp_norm=sqrt(Uexp(:,:,:,1).^2+Uexp(:,:,:,2).^2+Uexp(:,:,:,3).^2); 
Res=Ufit-Uexp;  % residual
Res_norm=Res(:,:,:,1).^2+Res(:,:,:,2).^2+Res(:,:,:,3).^2;

% % Plot the force field
% figure 
% simple_quiver(Ffit,'b');
% axis equal

% % Overlay the average field used for regularization:
% hold on
% simple_quiver(F_avg,'r');

% Plot the displacement fields
figure
hold off
simple_quiver(Uexp,'b')
hold on
simple_quiver(Ufit,'r')
axis equal

% % Plot the force field directly found from the experimental displacement:
% tmp=cat(4,cell_position,cell_position,cell_position);
% Uexp(tmp)=0;
% Fexp=M*Uexp(:);
% Fexp=reshape(Fexp,nx,ny,nz,3);
% %Pexp=dx*dy*dz*Fexp/(dy*dz); % dx*dy*dz*Ffit is the total force in a box. dxdy the amount of surface per box in um^2. The ratio is therefore in pN/um^2 i.e. Pascals. 
% figure
% simple_quiver(Fexp,'g')
% axis equal

% Plot the norm comparison
figure
plot(Uexp_norm(:),Ufit_norm(:),'*')
hold on
plot([0 max(Uexp_norm(:))],[0 max(Uexp_norm(:))])

% % % Plot the xz direction
% % figure
% % quiver(squeeze(Px(:,12,:)),squeeze(Pz(:,12,:)),squeeze(Uexp(:,12,:,1)),squeeze(Uexp(:,12,:,3)),'b')
% % hold on
% % quiver(squeeze(Px(:,12,:)),squeeze(Pz(:,12,:)),squeeze(Ufit(:,12,:,1)),squeeze(Ufit(:,12,:,3)),'r')
% % axis equal

% % Radius of the fitted sphere:
% mid=area_size/2;
% r2=(Px(cell_surface)*pxsize_xy-mid).^2+(Py(cell_surface)*pxsize_xy-mid).^2+(Pz(cell_surface)*pxsize_z-mid).^2;
% r=sqrt(r2);

% disp('Average force')
% disp(mean(Ffit_norm*pxsize_xy^2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Computation of the stress tensor sigma %%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=sparse(n_block,n_block);  % Used as a zero block
% Matrices that derive the components of the strain tensor epsilon (eps)
% from the displacement field U (with partial derivatives on the left l and right r). 
M_eps_xx_r=[Drx Z Z];
M_eps_xx_l=[Dlx Z Z];
M_eps_yy_r=[Z Dry Z];
M_eps_yy_l=[Z Dly Z];
M_eps_zz_r=[Z Z Drz];
M_eps_zz_l=[Z Z Dlz];
M_eps_xy_rr=1/2*[Dry Drx Z];
M_eps_xy_rl=1/2*[Dry Dlx Z];
M_eps_xy_lr=1/2*[Dly Drx Z];
M_eps_xy_ll=1/2*[Dly Dlx Z];
M_eps_xz_rr=1/2*[Drz Z Drx];
M_eps_xz_rl=1/2*[Drz Z Dlx];
M_eps_xz_lr=1/2*[Dlz Z Drx];
M_eps_xz_ll=1/2*[Dlz Z Dlx];
M_eps_yz_rr=1/2*[Z Drz Dry];
M_eps_yz_rl=1/2*[Z Drz Dly];
M_eps_yz_lr=1/2*[Z Dlz Dry];
M_eps_yz_ll=1/2*[Z Dlz Dly];
% Actual computation of the values of the strain tensor: 
% (The derivatives are computed on the right and left and averaged. 
% When the derivative is only available on one side, then this is
% kept as the value).
%Uexp(cell_position)=nan(sum(cell_position(:)),1);
eps_xx=nanmean([M_eps_xx_r*Ufit(:), M_eps_xx_l*Ufit(:)],2);
eps_yy=nanmean([M_eps_yy_r*Ufit(:), M_eps_yy_l*Ufit(:)],2);
eps_zz=nanmean([M_eps_zz_r*Ufit(:), M_eps_zz_l*Ufit(:)],2);
eps_xy=nanmean([M_eps_xy_rr*Ufit(:), M_eps_xy_rl*Ufit(:), M_eps_xy_lr*Ufit(:), M_eps_xy_ll*Ufit(:)],2);
eps_xz=nanmean([M_eps_xz_rr*Ufit(:), M_eps_xz_rl*Ufit(:), M_eps_xz_lr*Ufit(:), M_eps_xz_ll*Ufit(:)],2);
eps_yz=nanmean([M_eps_yz_rr*Ufit(:), M_eps_yz_rl*Ufit(:), M_eps_yz_lr*Ufit(:), M_eps_yz_ll*Ufit(:)],2);
eps_xx=reshape(eps_xx,nx,ny,nz);
eps_yy=reshape(eps_yy,nx,ny,nz);
eps_zz=reshape(eps_zz,nx,ny,nz);
eps_xy=reshape(eps_xy,nx,ny,nz);
eps_xz=reshape(eps_xz,nx,ny,nz);
eps_yz=reshape(eps_yz,nx,ny,nz);

% Put all the values together to have the strain as a tensor field, and
% derive the stress tensor field from it. 
strain=zeros(nx,ny,nz,3,3); % The last two indices are the indices of the strain matrix (xyz x xyz)
stress=zeros(nx,ny,nz,3,3); % The last two indices are the indices of the stress matrix (xyz x xyz)
principal_stress=zeros(nx,ny,nz);
principal_direction=zeros(nx,ny,nz,3);
strain(:,:,:,1,1)=eps_xx;
strain(:,:,:,1,2)=eps_xy;
strain(:,:,:,1,3)=eps_xz;
strain(:,:,:,2,1)=eps_xy;
strain(:,:,:,2,2)=eps_yy;
strain(:,:,:,2,3)=eps_yz;
strain(:,:,:,3,1)=eps_xz;
strain(:,:,:,3,2)=eps_yz;
strain(:,:,:,3,3)=eps_zz;
for i=1:nx
    for j=1:ny
        for k=1:nz
                strain_tmp=squeeze(strain(i,j,k,:,:));
                stress_tmp=2*mu*strain_tmp+lambda*trace(strain_tmp)*eye(3);
                stress(i,j,k,:,:)=stress_tmp;
                if ~isnan(sum(stress_tmp(:)))
                    [eigenvectors,eigenvalues]=eig(stress_tmp);
                    [v,p]=max(abs(diag(eigenvalues)));
                    principal_stress(i,j,k)=v;
                    principal_direction(i,j,k,:)=eigenvectors(:,p);
                else
                    principal_stress(i,j,k)=nan;
                end
        end
    end
end

% Compute the value of the pressure (stress*normal) on the surface of the cell: 
if model == 1
    stress_tmp=reshape(stress,nx*ny*nz,3,3);
    stress_cell_surface=stress_tmp(cell_surface(:),:,:);
    P=zeros(size(stress_cell_surface,1),1);
    if direction == 1
        norm_tmp=[1; 0; 0];
    elseif direction == 2
        norm_tmp=[0; 1; 0];
    elseif direction == 2
        norm_tmp=[0; 0; 1];
    end
    for i=1:size(stress_cell_surface,1)
        P(i)=-norm_tmp'*squeeze(stress_cell_surface(i,:,:))*norm_tmp;
    end
end

if model == 2
    stress_tmp=reshape(stress,nx*ny*nz,3,3);
    stress_cell_surface=stress_tmp(cell_surface(:),:,:);
    P=zeros(size(stress_cell_surface,1),1);
    if direction == 1
        norm_tmp=[0; 1; 0];
    elseif direction == 2
        norm_tmp=[0; 0; 1];
    elseif direction == 1
        norm_tmp=[0; 0; 1];
    end
    for i=1:size(stress_cell_surface,1)
        P(i)=-norm_tmp'*squeeze(stress_cell_surface(i,:,:))*norm_tmp;
    end
end

if model == 3
    % Define normal vectors (just radial unit vectors for the sphere):
    norm_x=(Px-n_px(1)/2);
    norm_y=(Py-n_px(2)/2);
    norm_z=(Pz-n_layers/2);
    norm_x=norm_x(cell_surface);
    norm_y=norm_y(cell_surface);
    norm_z=norm_z(cell_surface);
    tmp=sqrt(norm_x.^2+norm_y.^2+norm_z.^2);
    norm_x=norm_x./tmp;
    norm_y=norm_y./tmp;
    norm_z=norm_z./tmp;
    % To check normal vectors:
    %quiver3(Px(cell_surface(:)),Py(cell_surface(:)),Pz(cell_surface(:)),norm_x,norm_y,norm_z)
    
    % Compute normal pressure from stress tensor and normal vectors:
    stress_tmp=reshape(stress,nx*ny*nz,3,3);
    stress_cell_surface=stress_tmp(cell_surface(:),:,:);
    P=zeros(size(stress_cell_surface,1),1);
    for i=1:size(stress_cell_surface,1)
        norm_tmp=[norm_x(i); norm_y(i); norm_z(i)];
        P(i)=-norm_tmp'*squeeze(stress_cell_surface(i,:,:))*norm_tmp;
    end
    
    % Compute how aligned is the stress with the normal vector (theoretically should be aligned, but might be lost with noisy data or poor meshing, and can be recovered at least partially with the smoothing parameter):
    alignment=zeros(size(stress_cell_surface,1),1);
    stress_vector=zeros(3,size(stress_cell_surface,1));
    stress_vector_direction=zeros(3,size(stress_cell_surface,1));
    for i=1:size(stress_cell_surface,1)
        norm_tmp=[norm_x(i); norm_y(i); norm_z(i)];
        stress_vector(:,i)=-squeeze(stress_cell_surface(i,:,:))*norm_tmp; % Stress vector (Pa)
        stress_vector_direction(:,i)=stress_vector(:,i)/norm(stress_vector(:,i));    % Normalized stress vector (norm 1)
        alignment(i)=stress_vector_direction(:,i)'*norm_tmp;   % Value from -1 (opposite direction), 0 (orthogonal) to 1 (perfectly aligned). 
    end
end

% Plot the stress field
figure 
quiver3(Px(cell_surface),Py(cell_surface),Pz(cell_surface),stress_vector(1,:)',stress_vector(2,:)',stress_vector(3,:)','b')
axis equal
% hold on
% quiver3(Px(cell_surface),Py(cell_surface),Pz(cell_surface),stress_vector_direction(1,:)',stress_vector_direction(2,:)',stress_vector_direction(3,:)','g')

disp(' ')
toc
disp(' ')

disp('Ufit_max  Uexp_max  mean_%_error')
disp([max(Ufit_norm(:)) max(Uexp_norm(:)) 100*mean(abs((Uexp_norm(:)-Ufit_norm(:))./Uexp_norm(:)),'omitnan')])

if model==3
    disp('Average alignment of the stress vectors with the normal vectors and std')
    disp([mean(alignment(:),'omitnan') std(alignment(:),'omitnan')])
end

disp('Average pressure on the surface')
disp(mean(P))
disp('Std of the pressure on the surface')
disp(std(P))
