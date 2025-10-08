%% Visualisation VISJON-LP

clear all%, clf

mkdir FIGURES

cd COLORMAPS
    mat = dir('*.mat');
    for i=1:length(mat)
        load(mat(i).name);
    end    
cd ..


for num = [2]
    
    cd OUTPUT
    if exist(['Shearing_',num2str(1e6+num),'.mat'])
        load (['Shearing_',num2str(1e6+num),'.mat']);
        cd ..
    else
        cd ..
        break;
    end
    
    
    figure(11), clf, subplot(211)
    
        for i=1:20
            Phase1  =   find(Im == i);
            marksize = 1;
            hold on
            plot(xm(Phase1),ym(Phase1),'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end

    shading interp
    axis image, axis ij%, axis([100e3 400e3 0 50e3])
    title(['Time: ',num2str((time)/SecYear/1e6),' Ma'])
%     plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%     quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
    hold on
%     comp=['Himalaya_',num2str(num+1e6)];
%     print('-r150', '-djpeg', comp)
    ind=find(rho_s<2000);
    eta_s(ind)=NaN;
    E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
    eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
    T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
%     rho_s_2d    =   reshape(rho_s,ny+1,nx+1);
    strain_2d   =   reshape(strain,ny+1,nx+1);
%     D_2d        =   reshape(D_grain,ny+1,nx+1);
    
    xlabel('width [m]')
    ylabel('depth [m]')

    
    
    hold off
    
    subplot(212)
%     figure(2), clf
    colormap(flip(hawaii,1))
    
%     ind = find(rho_s_2d < 200);
%     E2nd_s_2d(ind)  = 1e-14;   
%     T2nd_s_2d(ind)  = 1e20;   
%     eta_s_2d(ind)   = 1e21;   
    
%         subplot(311)
    pcolor(x,y,log10(eta_s_2d(1:end-1,1:end-1))), caxis([17 25])
    shading interp
    colorbar
    axis image, axis ij%, axis([100e3 400e3 0 50e3])
    title('\eta [Pa.s]')
    hold on
%     contour(x,y,(log10(D_2d(1:end-1,1:end-1).*1000)),[-2:2],'k')
%     fill([0 0 1000e3 1000e3 0],[0 10e3 10e3 0 0],'w','EdgeColor','w'), hold on
%     plot([0 0 1000e3 1000e3 0],[0 670e3 670e3 0 0],'k-')
    xlabel('width [m]')
    ylabel('depth [m]')

%     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)

%     
%     subplot(223)
%     pcolor(x,y,T2nd_s_2d(1:end-1,1:end-1))%, caxis([18 24])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title('\tau [Pa]')
%     hold on
% %     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)
% %    
% %     
%     subplot(224)
%     pcolor(x,y,log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-17 -12])
%     shading interp
%     colorbar
%     axis image, axis ij
%     title('\dotepsilon [1/s]')
%     hold on
% %     fill([0 x Lx 0],[0 surface_y(1:5:end) 0 0],'w','EdgeColor','w','FaceColor','w','LineWidth',0.001)
% % 
%     
cd FIGURES
comp=['CompVisc_',num2str(num+1e6)];
print('-r300', '-djpeg', comp)
% %
cd ..

end