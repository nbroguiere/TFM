
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the contraction center maximizing the contraction and comparing with the centroid position:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Finding the contraction center maximizing the contraction and comparing with the centroid position');
O=zeros(3,n_timepoints);
O_center=zeros(3,n_timepoints);
contraction=zeros(n_timepoints,1);
O_center_index=zeros(3,n_timepoints);
contraction_at_centroid=zeros(n_timepoints,1);
for t=timepoints_to_process
    % Centroid position:
    bw=img_threshold(:,:,:,t);
    stats=regionprops('table',bw,'centroid','area');
    [value,pos]=max(stats.Area);
    O_center(:,t)=stats.Centroid(pos,:);
    
    % Contraction maximum:
    indices=cell_surface(:,:,:,t);
    
    tmp=Ffit(:,:,:,t,1);
    Fxt=tmp(indices(:));
    tmp=Ffit(:,:,:,t,2);
    Fyt=tmp(indices(:));
    tmp=Ffit(:,:,:,t,3);
    Fzt=tmp(indices(:));

    tmp=Px(:,:,:,t);
    Pxt=tmp(indices(:));
    tmp=Py(:,:,:,t);
    Pyt=tmp(indices(:));
    tmp=Pz(:,:,:,t);
    Pzt=tmp(indices(:));
       
    % Compute the total contraction force on the grid:
    contraction_matrix=zeros(nx,ny,nz);
    for i=1:nx
        for j=1:ny
            for k=1:nz
                Ox=Px(i,j,k,t);
                Oy=Py(i,j,k,t);
                Oz=Pz(i,j,k,t);
                contraction_matrix(i,j,k)=-sum((((Pxt-Ox).*Fxt)+(Pyt-Oy).*Fyt+(Pzt-Oz).*Fzt)./sqrt((Pxt-Ox).^2+(Pyt-Oy).^2+(Pzt-Oz).^2));
                contraction_matrix_xy(i,j,k)=-sum((((Pxt-Ox).*Fxt)+(Pyt-Oy).*Fyt)./sqrt((Pxt-Ox).^2+(Pyt-Oy).^2+(Pzt-Oz).^2));
                %                 contraction_matrix(i,j,k)=-sum((((Pxt-Ox).*Fxt)+(Pyt-Oy).*Fyt+(Pzt-Oz).*Fzt)./sqrt((Pxt-Ox).^2+(Pyt-Oy).^2+(Pzt-Oz).^2));
            end
        end
    end
    
    % Compute the contraction value at the centroid position
    if O_center(1,t)<xrange(1,end)
        O_center_index(1,t)=find(xrange>O_center(1,t),1);
    else
        O_center_index(1,t)=length(xrange);
    end
    if O_center(2,t)<yrange(1,end)
        O_center_index(2,t)=find(yrange>O_center(2,t),1);
    else
        O_center_index(2,t)=length(yrange);
    end
    if O_center(3,t)<zrange(1,end)
        O_center_index(3,t)=find(zrange>O_center(3,t),1);
    else
        O_center_index(3,t)=length(zrange);
    end
    contraction_at_centroid(t)=contraction_matrix(O_center_index(1,t),O_center_index(2,t),O_center_index(3,t));
    
    [contraction(t), pos]=max(contraction_matrix(:));
    [contraction_xy(t), posxy]=max(contraction_matrix_xy(:));
    O(1,t)=Px(pos);
    O(2,t)=Py(pos);
    O(3,t)=Pz(pos);
end
clear A b Pxt Pyt Pzt Fxt Fyt Fzt contraction_matrix

subplot(3,2,[1,3,5])
hold off
plot3(O(1,:),O(2,:),O(3,:),'b');
hold on
plot3(O(1,:),O(2,:),O(3,:),'*','color','b');
plot3(O_center(2,:),O_center(1,:),O_center(3,:),'r');
plot3(O_center(2,:),O_center(1,:),O_center(3,:),'*','color','r');
axis ij
%axis([0 1000 0 1000]);
title('Position mass vs contraction center');

subplot(3,2,2)
hold off
plot(O(1,:),'b');
hold on
plot(O_center(2,:),'r');
title('x position mass vs contraction center');

subplot(3,2,4)
V=zeros(3,n_timepoints);
V(:,2:end)=O_center(:,2:end)-O_center(:,1:end-1);
Vn=zeros(1,n_timepoints);
Vn(:)=sqrt(V(1,:).^2+V(2,:).^2);%+V(3,:).^2;
plot(Vn(1,:),'g');
title('Speed xy');

subplot(3,2,6)
plot(contraction_xy)
title('contraction xy');

set(gcf,'visible','on');
set(gcf,'Color',[1 1 1]);
if ~exist('output','dir')
    mkdir('output')
end
savefig('output/CentroidAndContraction.fig');
close(gcf);

figure
plot(contraction_at_centroid(timepoints_to_process))
title('Contraction at the centroid position');
set(gcf,'visible','on');
set(gcf,'Color',[1 1 1]);
savefig('output/Contraction_at_centroid.fig');
close(gcf);

clear contraction_matrix bw
toc
