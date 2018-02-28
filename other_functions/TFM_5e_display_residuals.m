% For which timepoint do you want to check the residuals:
t=30;

Ufit_norm_t=zeros(nx,ny,nz);
UexpC_norm_t=zeros(nx,ny,nz);
residual=zeros(nx,ny,nz);

for i=1:nx
    for j=1:ny
        for k=1:nz
            Ufit_norm_t(i,j,k)=norm(squeeze(Ufit(i,j,k,t,1:2)));
            UexpC_norm_t(i,j,k)=norm(squeeze(UexpC(i,j,k,t,1:2)));
            residual(i,j,k)=norm(squeeze(Ufit(i,j,k,t,1:2)-UexpC(i,j,k,t,1:2)));
        end
    end
end

figure
subplot(1,2,1)
plot(UexpC_norm_t(:),Ufit_norm_t(:),'.');
axis([0,3,0,3])
xlabel('|Uexp|')
ylabel('|Ufit|')

subplot(1,2,2)
plot(residual(:),'.');
xlabel('datapoints')
ylabel('residual')

if ~exist('output','dir')
    mkdir('output')
end

savefig('output/Residuals.fig');
