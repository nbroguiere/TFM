% The derived parameters need to be completely reconsidered with the Units.

% Need to redefine the cell surface as the elements containing the border.
% Better this way, otherwise I get a surface entirely bigger than the cell.

dx=step_size_xy*pxsize_xy;  % in um
dy=step_size_xy*pxsize_xy;  % in um
dz=step_size_z*pxsize_z;    % in um

% Calculation of the Lamé coefficients:
lambda=2*G*nu/(1-2*nu);   % lambda is related to the shear modulus and poisson ratio
mu=G;                     % mu identifies with the shear modulus G. 

% Adimensionalization of the forces through their weight parameters: 
F0=mu/(dx*dy*dz)^(1/3);                                   % Typical volumic force used for adimensionalization.
alpha_F_cell_surfaceAd=alpha_F_cell_surface/F0^2;         % Relative importance of minimizing the forces on the cell surface 
alpha_F_matrixAd=alpha_F_matrix/F0^2;                     % Relative importance of minimizing the forces compared to fitting the data in the matrix. Important note: in trials, when going up to 1e14, all fine. When going to 1e15, very different solution. Probably problems with the roundings. 
alpha_smoothingAd=alpha_smoothing/F0^2;                   % Relative importance of minimizing on the cell surface the difference between the force on a position and the average force on neighboring positions.

% Create the arrays storing fitted displacement and corresponding forces if not already existant. Units: px and pascals. 
disp('Finding the displacement field generated from forces restricted to the cell surface that matches the experimental displacement field best');
if ~exist('Ufit','var')
    Ufit=zeros(size(UexpC));
end
if ~exist('Ffit','var')
    Ffit=zeros(size(UexpC));
end

tic
for t=timepoints_to_process
    disp(['t=' num2str(t)]);
    disp(   'Building the matrix M of the mechanical model so that Force=M*Displacement');
    %%%%%%%%% Build the matrix giving the force from the displacement (in vector form):
    % Generate template matrices to build the blocks: 
    % (always put 0 in the border conditions: equivalent to zero padding the
    % data by zeros, reasonable since forces go to 0 far from the cell):
    n_block=nx*ny*nz;
    Identity=spdiags(ones(n_block,1),0,n_block,n_block);
    next=1;
    prev=-1;
    e=ones(n_block-1,1);
    e(nx:nx:n_block)=0;
    Mxn=spdiags([0; e],next,n_block,n_block);   % Mxn is the matrix that picks up the next element in the x direction
    Mxp=spdiags([e; 0],prev,n_block,n_block);
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
    
    % Fix the border conditions on the surface of the cells to also have reflected values: 
    indices=find(cell_surface(:,:,:,t));
    
    %cell_incl_surface=or(cell_position,cell_surface);
    for index=indices'
        [i,j,k]=ind2sub([nx,ny,nz],index);
        % Find the direction of the surface in x and fix what needs to be:
        if i==13 && j==15 && k==6 && t==30
            disp('debug');
        end
        if i==nx
        elseif cell_position(i+1,j,k,t)   % If the cell is just next, calls for the next value should return the local value. 
            Mxn(i+nx*(j-1)+nx*ny*(k-1),i+1+nx*(j-1)+nx*ny*(k-1))=0;
            Mxn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
        if i==1
        elseif cell_position(i-1,j,k,t)   % If the cell is just before, calls for the previous value should return the local value. 
            Mxp(i+nx*(j-1)+nx*ny*(k-1),i-1+nx*(j-1)+nx*ny*(k-1))=0;
            Mxp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
        % Same in y:
        if j==ny
        elseif cell_position(i,j+1,k,t)   % If the cell is just next, calls for the next value should return the local value. 
            Myn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j)+nx*ny*(k-1))=0;
            Myn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
        if j==1
        elseif cell_position(i,j-1,k,t)   % If the cell is just before, calls for the previous value should return the local value. 
            Myp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-2)+nx*ny*(k-1))=0;
            Myp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
        % Same in z:
        if k==nz
        elseif cell_position(i,j,k+1,t)   % If the cell is just next, calls for the next value should return the local value. 
            Mzn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k))=0;
            Mzn(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
        if k==1
        elseif cell_position(i,j,k-1,t)   % If the cell is just before, calls for the previous value should return the local value. 
            Mzp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-2))=0;
            Mzp(i+nx*(j-1)+nx*ny*(k-1),i+nx*(j-1)+nx*ny*(k-1))=1;
        end
    end
    
    % Then build the convolution filters h corresponding to the various second derivatives:
    second_derivative_straight=[1 -2 1];
    second_derivative_crossed=1/2*[1 -1 0; -1 2 -1; 0 -1 1];
    h=cell(3,3);
    h{1,1}=zeros(3,3,3);
    h{1,1}(:,2,2)=second_derivative_straight/(dx*dx);   % in um^-2
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

    % Put them all together to build the nine blocks to assemble to get the final matrix M satisfying F=M*U:
    % In this notation, i represents x y or z. 
    % block{i,j} gives the component of Fi coming from Uj.
    % Overall, the physical equation is 
    %     eps=1/2*(du_i/dx_j+du_j/dx_i)
    %     sigma=2*mu*eps+lambda*tr(eps)*I
    %     F=-sum(dsigma_ij/dx_i)
    block=cell(3,3);    % lambda and mu are in pN/um^2 and H in um^-2 so M is in pN/um^4
    for i=1:3
        for j=1:3
            if i==j
                block{i,j}=-(lambda+2*mu)*H{i,i};   % in pN/um^4
                for k=1:3
                    if k~=i
                        block{i,j}=block{i,j}-mu*H{k,k};    % in pN/um^4
                    end
                end
            else
                block{i,j}=-(lambda+mu)*H{i,j};   % in pN/um^4
            end
        end
    end
    M=[block{1,1} block{1,2} block{1,3}; block{2,1} block{2,2} block{2,3}; block{3,1} block{3,2} block{3,3}]; % in pN/um^4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/ building force matrix

    %%%%%%%% Similarly, build the average-of-neighboring-forces matrix %%%%%%%%
    % Define the standard convolution filter defining neighbours (here: 27),
    % not normalized (make it later on in the circle to take into account the
    % actual number of neighbours for each surface position line). 
    disp(   'Building the smoothing matrix');
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
    
    disp(   'Converting the experimental data to vector form');
    % Create V which is the field of displacement U, but in 1D vector-storage form:
    Vexp=UexpC(:,:,:,t,:);  % in um
    Vexp=Vexp(:);
    n_U=nx*ny*nz*3;
    
    disp(   'Creating the weight matrices');
    % Create a weight matrix with zeroes when the cross-correlation
    % coefficient has been deemed unacceptable:
    acceptable=~isnan(Vexp);
    % Then remove the nans from the experimental data, the data wont be used but nans would propagate even when multiplied by 0:
    Vexp(isnan(Vexp))=0;
    % (also apply a coefficient for how important is good fit of z data compared
    % to x and y data (note that z components of vectors are not as good
    % because of the quality of the raw data is far lesser):
    tmp=double(acceptable);
    tmp(n_block*2+1:n_block*3)=tmp(n_block*2+1:n_block*3)*sqrt(alpha_z);
    A=spdiags(tmp,0,sparse(n_U,n_U));
    
    % Create the coefficient matrix B for the force minimization (min of F'*B'*B*F=U'*M'*B'*B*M*U=norm(f(U)))
    % Note that B includes different coefficients for the matrix (alpha high, forces forbidden) and on the cell surface (alpha low, forces authorized)
    % For cell surface and borders authorize forces (excluding borders which are inside the cell):
    surface_and_borders=or(cell_surface(:,:,:,t),and(borders,~cell_position(:,:,:,t)));   % cell surface and borders which are not in the cell.
    if alpha_F_cell_surfaceAd>=0
        tmp=surface_and_borders(:,:,:)*sqrt(alpha_F_cell_surfaceAd);
    else
        tmp=surface_and_borders(:,:,:)*sqrt(-alpha_F_cell_surfaceAd);
    end    
    tmp2=[tmp(:);tmp(:);tmp(:)];
    B1=spdiags(tmp2(:),0,sparse(n_U,n_U));
    % For matrix only:
    tmp=~surface_and_borders(:,:,:)*sqrt(alpha_F_matrixAd);
    tmp2=[tmp(:);tmp(:);tmp(:)];
    B2=spdiags(tmp2(:),0,sparse(n_U,n_U));
    B=B1+B2;
    
    % Create the smoothing matrix D for the force regularization on the cell surface (min of F'*C'*C*F=U'*M'*C'*C*M*U=norm(F-Favg))
    tmp=cell_surface(:,:,:,t);              % Pick up the cell surface positions at this timepoint.
    tmp=sparse(tmp(:)');                    % Linearize the image to a vector
    H_avg_norm=tmp'.*H_avg.*tmp;            % Take a contribution coefficient of 0 for the positions not on the cell surface, and dont send values to positions outside the cell surface.
    H_avg_norm=H_avg_norm./max(sum(H_avg_norm,2),eps);  % Normalize each line so the average is actually an average (e.g. if there are 3 neighours, divide by 3). When there is a division by zero, divide by epsilon instead. 
    tmp=Identity-H_avg_norm;
    D=sqrt(alpha_smoothingAd)*blkdiag(tmp,tmp,tmp);
    
    disp(   'Solving the linear system');
    % Construct the matrices in a way that makes sure matlab recognizes
    % them as Hermitian for the linear solving: 
    tmp=A*Vexp;
    tmp2_1=B1*M;
    tmp2_2=B2*M;
    tmp3=D*M;
    if alpha_F_cell_surfaceAd>=0
        tmp4=A+tmp2_1'*tmp2_1+tmp2_2'*tmp2_2+tmp3'*tmp3;
    else
        tmp4=A-tmp2_1'*tmp2_1+tmp2_2'*tmp2_2+tmp3'*tmp3;      % This way of constructing the matrix makes sure matlab recognizes it as Hermitian, makes solving the system faster. 
    end
    Vfit=tmp4\tmp;              % Solving the finite element model
    Ffit_vector=M*Vfit;         % Computing the associated force field in pN/um^3
    Ufit(:,:,:,t,:)=reshape(Vfit,nx,ny,nz,3);  % Puting back from vector form to normal 3D vector field organization
    Ffit(:,:,:,t,:)=reshape(Ffit_vector,nx,ny,nz,3);
    toc
end
Ufitx=Ufit(:,:,:,:,1);
Ufity=Ufit(:,:,:,:,2);
Ufitz=Ufit(:,:,:,:,3);
Ufitx(cell_position)=nan(size(Ufitx(cell_position)));
Ufity(cell_position)=nan(size(Ufity(cell_position)));
Ufitz(cell_position)=nan(size(Ufitz(cell_position)));
Ufit=cat(5,Ufitx,Ufity,Ufitz);
clear A tmpx tmpy tmpz Vexp Ufit_solved Ffit_solved tmp tmp2_1 tmp2_2 tmp3 tmp4 A sqrtB V Ufit_tmp Ufitx Ufity Ufitz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computation of the strain and stress tensors epsilon and sigma %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Matrices to compute the first derivatives:
% Redefine the next and previous matrices with no cell edge specific border conditions:
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
I=Identity;
% derivative on the right: 
Drx=(Mxn-I)/dx;
Dry=(Myn-I)/dy;
Drz=(Mzn-I)/dz;
% derivative on the left:
Dlx=(I-Mxp)/dx;
Dly=(I-Myp)/dy;
Dlz=(I-Mzp)/dz;

% Matrices that derive the components of the strain tensor epsilon (eps)
% from the displacement field U (with partial derivatives on the left l and right r). 
Z=sparse(n_block,n_block);  % Used as a zero block
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
strain=zeros(nx,ny,nz,n_timepoints,3,3); % The last two indices are the indices of the strain matrix (xyz x xyz)
for t=timepoints_to_process
    Utmp=Ufit(:,:,:,t,:);
    strain(:,:,:,t,1,1)=reshape(nanmean([M_eps_xx_r*Utmp(:), M_eps_xx_l*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,2,2)=reshape(nanmean([M_eps_yy_r*Utmp(:), M_eps_yy_l*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,3,3)=reshape(nanmean([M_eps_zz_r*Utmp(:), M_eps_zz_l*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,1,2)=reshape(nanmean([M_eps_xy_rr*Utmp(:), M_eps_xy_rl*Utmp(:), M_eps_xy_lr*Utmp(:), M_eps_xy_ll*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,1,3)=reshape(nanmean([M_eps_xz_rr*Utmp(:), M_eps_xz_rl*Utmp(:), M_eps_xz_lr*Utmp(:), M_eps_xz_ll*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,2,3)=reshape(nanmean([M_eps_yz_rr*Utmp(:), M_eps_yz_rl*Utmp(:), M_eps_yz_lr*Utmp(:), M_eps_yz_ll*Utmp(:)],2),nx,ny,nz);
    strain(:,:,:,t,2,1)=strain(:,:,:,t,1,2);
    strain(:,:,:,t,3,1)=strain(:,:,:,t,1,3);
    strain(:,:,:,t,3,2)=strain(:,:,:,t,2,3);
end
clear M_eps_xx_r M_eps_xx_l M_eps_yy_r M_eps_yy_l M_eps_zz_r M_eps_zz_l M_eps_xy_rr M_eps_xy_rl M_eps_xy_lr M_eps_xy_ll M_eps_xz_rr M_eps_xz_rl M_eps_xz_lr M_eps_xz_ll M_eps_yz_rr M_eps_yz_rl M_eps_yz_lr M_eps_yz_ll

% Put all the values together to have the strain as a tensor field, and
% derive the stress tensor field from it. 
%
% Compute also the max principal stress (max eigenvalue)
%
% strain(i,j,k,t,i2,j2) has no unit. 
% stress(i,j,k,t,i2,j2) in Pa. 
%
stress=zeros(nx,ny,nz,n_timepoints,3,3); % The last two indices are the indices of the stress matrix (xyz x xyz)
principal_stress=zeros(nx,ny,nz,n_timepoints);
principal_direction=zeros(nx,ny,nz,n_timepoints,3);
for t=timepoints_to_process
    for i=1:nx
        for j=1:ny
            for k=1:nz
                strain_tmp=squeeze(strain(i,j,k,t,:,:));
                stress_tmp=2*mu*strain_tmp+lambda*trace(strain_tmp)*eye(3);
                stress(i,j,k,t,:,:)=stress_tmp;
                if ~isnan(sum(stress_tmp(:)))
                    [eigenvectors,eigenvalues]=eig(stress_tmp);
                    [v,p]=max(abs(diag(eigenvalues)));
                    principal_stress(i,j,k,t)=v;
                    principal_direction(i,j,k,t,:)=eigenvectors(:,p);
                else
                    principal_stress(i,j,k,t)=nan;
                end
            end
        end
    end
end
clear eigenvalues strain_tmp stress_tmp
