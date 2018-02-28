% Derived values using the mesh that are kept for data analysis:
positions=cell(n_timepoints,1); % in um 
normals=cell(n_timepoints,1);   % no unit (normalized to 1)
position_area=cell(n_timepoints,1);
cell_surface_area=zeros(n_timepoints,1);    % in square um
pressure_vector=cell(n_timepoints,1);
norm_of_pressure=cell(n_timepoints,1);
normal_pressure=cell(n_timepoints,1);
max_eigenvalue=cell(n_timepoints,1);
average_pressure=zeros(n_timepoints,1);
average_normal_pressure=zeros(n_timepoints,1);
average_norm_of_pressure=zeros(n_timepoints,1);
average_max_eigenvalue=zeros(n_timepoints,1);

% Create the output folders
if ~exist('output3D_tif','dir')
    mkdir('output3D_tif')
end
if ~exist('output3D_fig','dir')
    mkdir('output3D_fig')
end

key_frames(1)=1;
n_frames=key_frames(end);

i=1:n_frames;
time=interp1q(key_frames',key_times',i');
azimuth=interp1q(key_frames',key_azimuths',i');
elevation=interp1q(key_frames',key_elevations',i');
% interplq can have some rounding errors like e.g. -3.55e-15 that mess up
% integer values when they occur, so when the user obviously wants integer 
% values, do a round up: 
for i=1:n_frames
    if time(i)-round(time(i))<0.01
        time(i)=round(time(i));
    end
end

% Note that no check that these values are within the bounds! Maybe add it for more safety later on. 
% Nice view position: 46 52

Ufitx=Ufit(:,:,:,:,1);
Ufity=Ufit(:,:,:,:,2);
Ufitz=Ufit(:,:,:,:,3);

disp('3D rendering')
%screen_size = get(0,'ScreenSize'); 
screen_size = [100 100 1700 900]; 
indices_rand_multiplier=rand(nx*ny*nz,1)<plot_fraction;

% Pick up the thresholded image and smooth it:
if ~exist('img_threshold_smooth','var') || sigma_currently_used~=sigma
    disp('Smoothing the segmented cell surface')
    tic
    sigma_currently_used=sigma;
    img_threshold_smooth=zeros(size(img_threshold)); % Left with x y z t entries
    for t=1:n_timepoints
        % Pick up the stack for the timepoint under consideration:
        img_tmp=double(img_threshold(:,:,:,t));  
        % Smooth it (3D gaussian):
        box_range=-floor(3*sigma/pxsize_xy):floor(3*sigma/pxsize_xy);
        filter=normpdf(box_range,0,sigma/pxsize_xy);
        img_tmp=imfilter(img_tmp,filter);
        img_tmp=imfilter(img_tmp,filter');
        box_range=-floor(3*sigma/pxsize_z):floor(3*sigma/pxsize_z);
        filter=zeros(1,1,length(box_range));
        filter(:)=normpdf(box_range,0,sigma/pxsize_z);
        img_tmp=imfilter(img_tmp,filter);
        img_threshold_smooth(:,:,:,t)=img_tmp;
    end
    toc
    clear img_tmp
end

%Prepare everything for the colouring of the faces:
%force=Ffit(:,:,:,:,1).^2+Ffit(:,:,:,:,2).^2+Ffit(:,:,:,:,3).^2;   % Norm of the force on each point
%force=force./max(force);
for i=1:n_frames
    tic
    if i>1 && time(i)==time(i-1)    % No need to recompute anything if just a change of view.
        view(azimuth(i),elevation(i))
    else
        disp(['Frame=' num2str(i) '/' num2str(n_frames)])
        close(gcf);
        figure('Color',[0 0 0]);
        %set(gcf,'visible','off');
        fig = gcf;
        set(fig,'position',screen_size);
        set(fig,'InvertHardcopy','off');
        
        % Make a linear interpolation between the images at the two closest
        % timepoints (this is OK as soon as the images are close enough so that
        % the smoothing creates a smooth transition. Otherwise, it would need a
        % morphing instead).
        t=floor(time(i));
        dt=time(i)-floor(time(i));
        if dt==0
            img_tmp=permute(img_threshold_smooth(:,:,:,t),[2,1,3]); % need to permute i and j to keep the patch in the same coordinate system as the rest given matlab patch yx convention
        else
            img_tmp=permute((1-dt)*img_threshold_smooth(:,:,:,t)+dt*img_threshold_smooth(:,:,:,t+1),[2,1,3]);
        end
        
        % Find the surface
        [xx, yy, zz]=meshgrid((1:n_px(1))*pxsize_xy,(1:n_px(2))*pxsize_xy,(1:n_layers)*pxsize_z);
        fv=isosurface(xx,yy,zz,img_tmp,0.5);
        fv=reducepatch(fv,patch_reduction);
        n_vertices=size(fv.vertices,1);
        
        % Find the facet centers and normal unit vectors:
        a = fv.vertices(fv.faces(:,2),:)-fv.vertices(fv.faces(:,1),:);
        b = fv.vertices(fv.faces(:,3),:)-fv.vertices(fv.faces(:,1),:);
        c = cross(a, b, 2);
        area_faces=sqrt(sum(c.^2,2));
        normals_faces=-c./cat(2,area_faces,area_faces,area_faces);
        facet_centers=(fv.vertices(fv.faces(:,1),:)+fv.vertices(fv.faces(:,2),:)+fv.vertices(fv.faces(:,3),:))/3;   % coordinates in um. 
                % to check the normals, use:
                % quiver3(facet_centers(:,1),facet_centers(:,2),facet_centers(:,3),normals_faces(:,1),normals_faces(:,2),normals_faces(:,3),2)
                % quiver3(fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),normals_vertices(:,1),normals_vertices(:,2),normals_vertices(:,3),2)
                
         % Associate with each vertex the average of the normals of the connected faces:
        normals_vertices=zeros(n_vertices,3);
        associated_area=zeros(n_vertices,1);
        for m=i:n_vertices
            associated_faces=[find(fv.faces(:,1)==m); find(fv.faces(:,2)==m); find(fv.faces(:,3)==m)] ;
            normals_vertices(m,:)=mean(normals_faces(associated_faces,:),1);
            normals_vertices(m,:)=normals_vertices(m,:)/sqrt(normals_vertices(m,1)^2+normals_vertices(m,2)^2+normals_vertices(m,3)^2);
            associated_area(m)=sum(area_faces(associated_faces))/3;
        end
        % Save the position of the vertices and associated normals and surface areas:
        if ~(i>1 && floor(time(i))==floor(time(i-1)))  % Save only once for each physical timepoint, do not save for interpolated timepoints. 
            normals{t}=normals_vertices;
            positions{t}=fv.vertices;
            position_area{t}=associated_area;
        end
        
        % Positions and values of the strains known on the cell surface: 
        tmp=cell_surface(:,:,:,t);
        indices=find(tmp);
        
        tmp=squeeze(stress(:,:,:,t,:,:));
        tmp=reshape(tmp,nx*ny*nz,3,3);
        stress_t=tmp(indices,:,:);      % stress in Pa

        tmp=Px(:,:,:,t);
        Pxt=tmp(indices)*pxsize_xy;
        tmp=Py(:,:,:,t)*pxsize_xy;
        Pyt=tmp(indices);
        tmp=Pz(:,:,:,t)*pxsize_z;
        Pzt=tmp(indices);               % coordinates in um. 
        
        if dt>0
            tmp=cell_surface(:,:,:,t+1);
            indices2=find(tmp);
            
            tmp=squeeze(stress(:,:,:,t+1,:,:));
            tmp=reshape(tmp,nx*ny*nz,3,3);
            stress_t2=tmp(indices2,:,:);      % stress in Pa
            
            tmp=Px(:,:,:,t+1);
            Pxt2=tmp(indices2)*pxsize_xy;
            tmp=Py(:,:,:,t+1)*pxsize_xy;
            Pyt2=tmp(indices2);
            tmp=Pz(:,:,:,t+1)*pxsize_z;
            Pzt2=tmp(indices2);               % coordinates in um. 
        end
        
        % Remove from the known stress tensors list the few values that have a nan inside:
        tmp=sum(sum(stress_t,2),3);
        not_nan_positions=find(~isnan(tmp));
        Pxt=Pxt(not_nan_positions);
        Pyt=Pyt(not_nan_positions);
        Pzt=Pzt(not_nan_positions);
        stress_t=stress_t(not_nan_positions,:,:);
        
        if dt>0
            tmp=sum(sum(stress_t2,2),3);
            not_nan_positions=find(~isnan(tmp));
            Pxt2=Pxt2(not_nan_positions);
            Pyt2=Pyt2(not_nan_positions);
            Pzt2=Pzt2(not_nan_positions);
            stress_t2=stress_t2(not_nan_positions,:,:);
        end
        clear not_nan_positions
        
        % Gaussian smoothing and interpolation of the stress tensor on the cell surface:
        stress_t_smooth=zeros(n_vertices,3,3);   % smoothed stress tensors on the facet positions
        for m=1:n_vertices
            fc=fv.vertices(m,:);                                               % facet center considered, position is in um
            r2=(Pxt(:)-fc(1)).^2+(Pyt(:)-fc(2)).^2+(Pzt(:)-fc(3)).^2;          % square distance to known stresses
            nf=sum(exp(-r2/gaussian_smoothing_stress_display^2),'omitnan');    % normalization factor for the gaussian averaging
            stress_t_weighted=zeros(size(stress_t));
            for p=1:size(stress_t,1)
                stress_t_weighted(p,:,:)=stress_t(p,:,:)*exp(-r2(p)/gaussian_smoothing_stress_display^2);   % gaussian_smoothing_stress_display in um
            end
            stress_t_smooth(m,:,:)=sum(stress_t_weighted,1,'omitnan')/nf;                          % Gaussian averaged stress calculation
            if dt>0
                r2=(Pxt2(:)-fc(1)).^2+(Pyt2(:)-fc(2)).^2+(Pzt2(:)-fc(3)).^2;
                nf=sum(exp(-r2/gaussian_smoothing_stress_display^2),'omitnan');
                stress_t_weighted2=zeros(size(stress_t2));
                for p=1:size(stress_t2,1)
                    stress_t_weighted2(p,:,:)=stress_t2(p,:,:)*exp(-r2(p)/gaussian_smoothing_stress_display^2);   % gaussian_smoothing_stress_display in um
                end
                stress_t_smooth(m,:,:)=(1-dt)*stress_t_smooth(m,:,:)+dt*sum(stress_t_weighted2,1,'omitnan')/nf;
            end
        end
        clear r2 fc nf stress_t_weighted stress_t_weighted2

        % Compute on the surface of the cell the norm of the pressure, the value of the normal pressure, and max |eigenvalue| of the stress tensor(used for the color plot and for the average pressure calculation):
        pressure=zeros(n_vertices,1);
        pressure_vector_tmp=zeros(n_vertices,3);
        norm_of_pressure_tmp=zeros(n_vertices,1);
        normal_pressure_tmp=zeros(n_vertices,1);
        max_eigenvalue_tmp=zeros(n_vertices,1);
        for m=1:n_vertices
            pressure_vector_tmp(m,:)=-squeeze(stress_t_smooth(m,:,:))*[normals_vertices(m,1); normals_vertices(m,2); normals_vertices(m,3)];    % minus sign to get the pressure of the organoid on the matrix rather than the matrix on the organoid. 
            norm_of_pressure_tmp(m)=norm(pressure_vector_tmp(m,:));
            normal_pressure_tmp(m)=[normals_vertices(m,1); normals_vertices(m,2); normals_vertices(m,3)]'*pressure_vector_tmp(m,:)';
            max_eigenvalue_tmp(m)=max(abs(eig(squeeze(stress_t_smooth(m,:,:)))));
        end
        if type_of_color_display==0
            pressure=norm_of_pressure_tmp;
            color_scale_label='Pressure (Pa)';
        elseif type_of_color_display==1
            pressure=normal_pressure_tmp;
            color_scale_label='Normal pressure (Pa)';
        elseif type_of_color_display==2
            pressure=max_eigenvalue_tmp;
            color_scale_label='Principal stress (Pa)';
        end
        % To check the pressure vectors use:
        % quiver3(fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3),pressure_vector_tmp(:,1),pressure_vector_tmp(:,2),pressure_vector_tmp(:,3),2)
        
        % Compute the average pressure (integral of the pressure over the surface of the cell) if it wasn't computed yet:
        % In some very rare cases, e.g. a vertex connected to faces whose normals average to 0, 'nan' can appear. Make sure this doesn't crash the whole computation. 
        if cell_surface_area(t)==0  % as a way to detect that this timepoint has not yet been processed. 
            % Save the whole field of values:
            pressure_vector{t}=pressure_vector_tmp;
            norm_of_pressure{t}=norm_of_pressure_tmp;
            normal_pressure{t}=normal_pressure_tmp;
            max_eigenvalue{t}=max_eigenvalue_tmp;
            % Save the average values:
            first_vertex_of_each_face=fv.faces(:,1);
            cell_surface_area(t)=sum(area_faces,'omitnan');
            average_normal_pressure(t)=sum(normal_pressure_tmp(first_vertex_of_each_face).*area_faces,'omitnan')/cell_surface_area(t);
            average_norm_of_pressure(t)=sum(norm_of_pressure_tmp(first_vertex_of_each_face).*area_faces,'omitnan')/cell_surface_area(t);
            average_max_eigenvalue(t)=sum(max_eigenvalue_tmp(first_vertex_of_each_face).*area_faces,'omitnan')/cell_surface_area(t);
        end
        
        % Create the colormap and color values:
        tmp=jet(342);
        tmp=tmp(44:299,:);
        colormap(tmp);
        if normalize_color
            pressure_normalized=(pressure-minimum_color)/(maximum_color-minimum_color);
        else
            minimum_color=min(pressure);
            maximum_color=max(pressure);
            pressure_normalized=(pressure-minimum_color)/(maximum_color-minimum_color);
        end
        if inverted_color_display
            pressure_normalized=1-pressure_normalized;
        end
        color=pressure_normalized*size(colormap,1);
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Compute local value of the total contraction for coloring:
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Compute the total contraction force on the vertices position:
%         contraction_matrix=zeros(size(fv.vertices,1),1);
%         for j=1:size(fv.vertices,1)
%             Ox=fv.vertices(j,2)/pxsize_xy;
%             Oy=fv.vertices(j,1)/pxsize_xy;
%             Oz=fv.vertices(j,3)/pxsize_z;
%             contraction_matrix(j)=-sum((((Pxt-Ox).*Fxt)+(Pyt-Oy).*Fyt+(Pzt-Oz).*Fzt)./sqrt((Pxt-Ox).^2+(Pyt-Oy).^2+(Pzt-Oz).^2));
%         end
%         
%          % Contraction t+1:
%         if dt>0
%             tmp=cell_surface(:,:,:,t+1);
%             indices=tmp(:);    
%             
%             tmp=FxFit(:,:,:,t+1);
%             Fxt=tmp(indices(:));
%             tmp=FyFit(:,:,:,t+1);
%             Fyt=tmp(indices(:));
%             tmp=FzFit(:,:,:,t+1);
%             Fzt=tmp(indices(:));
%             
%             tmp=Px(:,:,:,t+1);
%             Pxt=tmp(indices(:));
%             tmp=Py(:,:,:,t+1);
%             Pyt=tmp(indices(:));
%             tmp=Pz(:,:,:,t+1);
%             Pzt=tmp(indices(:));
%             
%             % Compute the total contraction force on the vertices positions:
%             contraction_matrix2=zeros(size(fv.vertices,1),1);
%             for j=1:size(fv.vertices,1)
%                 Ox=fv.vertices(j,2)/pxsize_xy;
%                 Oy=fv.vertices(j,1)/pxsize_xy;
%                 Oz=fv.vertices(j,3)/pxsize_z;
%                 contraction_matrix2(j)=-sum((((Pxt-Ox).*Fxt)+(Pyt-Oy).*Fyt+(Pzt-Oz).*Fzt)./sqrt((Pxt-Ox).^2+(Pyt-Oy).^2+(Pzt-Oz).^2));
%             end
%             contraction_matrix=(1-dt)*contraction_matrix+dt*contraction_matrix2;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if ~normalize_color
%             %color=griddata(Py_tmp*pxsize_xy,Px_tmp*pxsize_xy,Pz_tmp*pxsize_z,force_tmp,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),'linear');
%             color=contraction_matrix./max(contraction(:))*size(colormap,1);  % Normalize to the overall maximum contraction
%         else
%             color=contraction_matrix/contraction_norm*size(colormap,1);
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %             Plot the cell surface:             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p=patch(fv,'FaceVertexCData',color,'FaceColor','interp','EdgeColor','none','CDataMapping','direct');
        daspect([1,1,1])
        axis([0,n_px(1)*pxsize_xy,0,n_px(2)*pxsize_xy,0,n_layers*pxsize_z])
        view(3);
        lighting gouraud
        set(gca,'Color',[0 0 0])
        set(gca,'XColor',[1 1 1])
        set(gca,'YColor',[1 1 1])
        set(gca,'ZColor',[1 1 1])
        view(azimuth(i),elevation(i))
%         if ~exist('light_handle','var')
%             light_handle=camlight('left');
%         else
%             camlight(light_handle,'left');
%         end
        camlight
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         Plot the displacement field:           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if dt==0
            Uxtmp=Ufitx(:,:,:,t);
            Uytmp=Ufity(:,:,:,t);
            Uztmp=Ufitz(:,:,:,t);    
        else
            Uxtmp=(1-dt)*Ufitx(:,:,:,t)+dt*Ufitx(:,:,:,t+1);
            Uytmp=(1-dt)*Ufity(:,:,:,t)+dt*Ufity(:,:,:,t+1);
            Uztmp=(1-dt)*Ufitz(:,:,:,t)+dt*Ufitz(:,:,:,t+1);
        end
        
        % Keep only a subpart of the vectors and scale them: 
        indices=~isnan(Uxtmp);
        indices=logical(indices(:).*indices_rand_multiplier);
        
        Pxt=Px(:,:,:,t);
        Pxt=Pxt(indices);
        Pyt=Py(:,:,:,t);
        Pyt=Pyt(indices);
        Pzt=Pz(:,:,:,t);
        Pzt=Pzt(indices);
        Uxt=Uxtmp(indices);
        Uyt=Uytmp(indices);
        Uzt=Uztmp(indices);
        
        hold on
        quiver3(Pxt*pxsize_xy,Pyt*pxsize_xy,Pzt*pxsize_z,Uxt*f_U_3D,Uyt*f_U_3D,Uzt*f_U,0,'w');
        xlabel('x, micrometers')
        ylabel('y, micrometers')
        zlabel('z, micrometers')
        title(['Time: ' num2str(floor(time(i)*timestep)) timeunit],'color','w');
        
        if minimum_color==0
            c=colorbar('Ticks',[1,256],'TickLabels',{num2str(0),num2str(maximum_color)});
        elseif minimum_color>0
            c=colorbar('Ticks',[1,256],'TickLabels',{num2str(minimum_color),num2str(maximum_color)});
        elseif minimum_color<0
            position_of_zero=1-255*minimum_color/(maximum_color-minimum_color);
            c=colorbar('Ticks',[1,position_of_zero,256],'TickLabels',{num2str(minimum_color),num2str(0),num2str(maximum_color)});
        end
        c.Label.String=color_scale_label;
        c.Label.Color='w';
        c.Color='w';
    end
    
    % Save the figure:
    savefig(['output3D_fig\frame' num2str(i) '.fig'])
    saveas(gcf,['output3D_tif\display4D_i' num2str(i) '.tif']);
    %%%%print(['output3D_tif\display4D_i' num2str(i) '.tif'],'-dtiff','-r100');
    %F = getframe(fig) % not sure I would have control of the resolution, but at least I wouldnt have to export individual pictures... 
    toc
end
clear Pxt Pyt Pzt Uxt Uyt Uzt

disp('The average norm of the pressure and normal pressure were computed at each timepoint under average_normal_pressure and average_norm_of_pressure');

% % Set back the figure display to normal parameters (visible, white background). 
% set(gcf,'visible','on');
% set(gcf,'Color',[1 1 1]);
% close(gcf);
