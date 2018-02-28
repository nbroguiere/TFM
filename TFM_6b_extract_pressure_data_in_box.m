xmin=user_box(1);
xmax=user_box(2);
ymin=user_box(3);
ymax=user_box(4);
zmin=user_box(5);
zmax=user_box(6);

P_box=zeros(n_timepoints,1);    % Pa which are pN/um^2
F_box=zeros(n_timepoints,1);    % pN
for t=timepoints_to_process
    pos_tmp=positions{t};
    
    indices=find(and(and(and(and(and(pos_tmp(:,1)>xmin,pos_tmp(:,1)<xmax),pos_tmp(:,2)>ymin),pos_tmp(:,2)<ymax),pos_tmp(:,3)>zmin),pos_tmp(:,3)<zmax));
    
    if pressure_definition==0
        P_tmp=norm_of_pressure{t}(indices);
    elseif pressure_definition==1
        P_tmp=normal_pressure{t}(indices);
    elseif pressure_definition==2
        P_tmp=max_eigenvalue{t}(indices);
    end
    
    area_tmp=position_area{t}(indices);
    P_box(t)=sum(P_tmp(:).*area_tmp(:),'omitnan')/sum(area_tmp(:));

    Pvector_tmp=pressure_vector{t}(indices,:);
    F_box_vector=area_tmp(:)'*Pvector_tmp;
    F_box(t)=norm(F_box_vector);
    
    quiver3(pos_tmp(indices,1),pos_tmp(indices,2),pos_tmp(indices,3),pressure_vector{t}(indices,1),pressure_vector{t}(indices,2),pressure_vector{t}(indices,3));
end
clear pos_tmp xmin xmax ymin ymax zmin zmax area_tmp indices Pvector_tmp F_box_vector

figure
subplot(1,2,1)
plot(timepoints_to_process*timestep,P_box(timepoints_to_process));
xlabel(['Time (' timeunit ')']);
if pressure_definition==0
    ylabel(['Average pressure (Pa) with box limits [' num2str(user_box) ']']);
elseif pressure_definition==1
    ylabel('Average normal pressure (Pa)');
elseif pressure_definition==2
    ylabel('Average maximum eigenvalue of the stress tensor (Pa)');
end


subplot(1,2,2)
plot(timepoints_to_process*timestep,F_box(timepoints_to_process));
xlabel(['Time (' timeunit ')']);
ylabel(['Total force (pN) within box limits']);

if ~exist('output','dir')
    mkdir('output')
end
savefig('output/Pressure in custom box.fig');
disp('The values of the average pressure in the box are stored in P_box (in Pa)');
disp('and the values of the total force in the box are stored in F_box (in pN)');
