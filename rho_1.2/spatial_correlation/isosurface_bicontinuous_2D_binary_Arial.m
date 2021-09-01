clear
close all

n_particle = 6912;
N = n_particle;       % number of particles
n_conf = 1;

%% time
dt = 1e-2;
q = 7.25;
n_ave = 1;
n_inteval = 1;
n_block = n_ave*n_inteval;

ts = [0 1:9];
for n_2 = 1:7
    ts = [ts n_ave*(10^n_2)*[1:9]];
end
ts = ts*n_inteval;
ts = ts(ts<=1e7);

n_frame = length(ts);

t = ts;

iT = 0;
% ic = 1;

index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
index_time = [12 21 30 39];

Temperature = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
Temperature = [1.0];
%% start loop
for T_i = fliplr(Temperature) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    %     mean_MC_c = zeros(length(index_time),1);
    %     mean_GC_c = zeros(length(index_time),1);
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        %% time origin
        %%%%%%%%%%%%%%%% load rg.mat file %%%%%%%%%%%%%%%%
        filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
        load(filename_rg)
        
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        [p_0,coord_0] = calculate_shape(rg_all);
        for i1 = 1:size(rg_all,3)
            [V,D] = eig(rg_all(:,:,i1));
            [lambda,ind] = sort(diag(D));
            Vs_0(:,:,i1) = V(:,ind);
        end
        
        %%%%%%%%%%%%%%%% load trajectory file %%%%%%%%%%%%%%%%
        filename = ['../',num2str(T_i,'%.1f'),'/eq.',num2str(ic,'%.0f'),'.dump'];
        dump = readmatrix(filename,'FileType','text');
        
        index = find(isnan(dump(:,3))==0);
        
        id = dump(index,1);
        r = dump(index,3:5);
        l = dump(1:3,1:2);
        
        type = dump(index,2);
        
        %% PBC
        Lx = (l(1,2)-l(1,1));
        Ly = (l(2,2)-l(2,1));
        Lz = (l(3,2)-l(3,1));
        L = [Lx, Ly, Lz];     % box size
        pbc = [1, 1, 1];      % boundary conditions
        r = r-round(r./L).*L;
        
        itime = 0;
        %% loop over time
        for it = index_time
            figure
            disp(['t = ', num2str(t(it))])
            itime = itime+1;
            rg_all = rg_all_t_c(:,:,:,it,ic);
            [p_t,coord_t] = calculate_shape(rg_all);
            for i1 = 1:size(rg_all,3)
                [V,D] = eig(rg_all(:,:,i1));
                [lambda,ind] = sort(diag(D));
                Vs(:,:,i1) = V(:,ind);
            end
            
            %%%%%%%%%%%%%%%% calculate correlation %%%%%%%%%%%%%%%%
            Ct = (3*(dot(Vs_0(:,3,:),Vs(:,3,:))).^2-1)/2;
            Ct = Ct(:);
            C_th = ((1-(1/2/pi))^2*3-1)/2;
%             C_th = mean(Ct);
            Ct_center = Ct-C_th;
            
            p = [Ct_center];
            
            %         % scatter
            %         figure(1)
            %         scatter3(X,Y,Z,50,p(index_selection),'filled')
            %         axis equal
            
            %% contour
            index_t0 = ((1-1)*N+(1:N));
            index_t = ((it-1)*N+(1:N));
            rt0 = r(index_t0,:);
            rt = r(index_t,:);
            
            %%%%%%%%%%%%%%%% define region %%%%%%%%%%%%%%%%
            thickness = 16;
            center = 0;
            
            origin_X = 0;
            origin_Y = 0;
            origin_Z = 0;
            
            %         index_selection = find(abs(rt0(:,3)-center)<thickness/2);
            index_selection = 1:N;
            X = rt0(index_selection,1);
            Y = rt0(index_selection,2);
            Z = rt0(index_selection,3);
            
            % shadow
            Xs = [];
            Ys = [];
            Zs = [];
            ps = [];
            for ix = -1:1
                for iy = -1:1
                    for iz = -1:1
                        Xs = [Xs;X+Lx*ix];
                        Ys = [Ys;Y+Ly*iy];
                        Zs = [Zs;Z+Lz*iz];
                        ps = [ps;p];
                    end
                end
            end
            
            edge = 2.5;
            index_shadow = find(abs((Xs)+origin_X)<abs(l(1,2)-l(1,1))/2+edge&...
                abs((Ys)+origin_Y)<abs(l(2,2)-l(2,1))/2+edge);
            
            Xs = Xs(index_shadow);
            Ys = Ys(index_shadow);
            Zs = Zs(index_shadow);
            ps = ps(index_shadow,:);
            
            dl = Lx/200; % pixel size
            n_edge = 1;
            [XX, YY, ZZ] = meshgrid((-n_edge*dl:dl:Lx+n_edge*dl)-Lx/2,(-n_edge*dl:dl:Ly+n_edge*dl)-Ly/2,3.5);
            F = scatteredInterpolant(Xs, Ys, Zs, ps, 'natural');
            CC = F(XX, YY, ZZ);
            %             CC = griddata(Xs, Ys, Zs, ps, XX, YY, ZZ, 'natural');
            IC = CC>0;
            f_C = sum(IC,'all')/length(XX(:)); % fraction of contious phase
            
            %% isosurface
            hold on
            
            CC(CC+C_th>1) = 1;
            CC(CC+C_th<-0.5) = -0.5;
            CC_r = CC;
            CC_r(CC>0) = (CC(CC>0)-min(CC(CC>0)))/(max(CC(CC>0))-min(CC(CC>0)));
            CC_r(CC<0) = (CC(CC<0))/(-min(CC(CC<0)));
            
%             pcolor(XX, YY, CC_r);
                       
            color = gray(256);
            colormap(color(64:256,:));
            
%             colormap([0 0.18 0.28; 1 1 1])

%             
            % outline
            contourf(-XX, -YY, CC, [-100 0 100],'-k','LineWidth',2);
            
            % contour_dark
            contour(-XX, -YY, CC, [0:-1/3:-1]*abs(min(CC,[],'all')),'-','Color','w','LineWidth',1);
            % contour_light
            contour(-XX, -YY, CC, [0:1/3:1]*abs(max(CC,[],'all')),'-','Color','k','LineWidth',1);

%             % contour
%             contour(XX, YY, CC, [1:-1/3:-1]*(max(CC,[],'all')-min(CC,[],'all')),'-','Color','#808080','LineWidth',1);

% %             color = jet(128);
% %             color = color([1:40,89:128],:);
%             
%             color_Y = [1,1,0] + [0,0,1].*[0:0.05:1]';
%             color_B = [0.1,0.4,1].*[0:0.05:1]';
%             color = [color_B;color_Y];
%             
%             load('YlOrRd.mat')
%             color = rgba(:,1:3);
% 
% %             colormap(flipud(color))
%             caxis([-1 1])
            
            set(gca,'YDir','normal')
            
            %% figure setting
            box on
            axis equal
            %         shading interp
            
            %         xlim([0,l(1,2)-l(1,1)])
            %         ylim([0,l(2,2)-l(2,1)])
            %         zlim([0,l(3,2)-l(3,1)])
            
            xlim([-Lx/2 Lx/2])
            ylim([-Ly/2 Ly/2])
            zlim([-Lz/2 Lz/2])
            
            xticks(-20:5:20)
            yticks(-20:5:20)
            
%             xlabel('x','FontSize',24,'Interpreter','latex')
%             ylabel('y','FontSize',24,'Interpreter','latex')
            
            xlabel('x','FontSize',24)
            ylabel('y','FontSize',24)
            
            set(gca,'LineWidth',2)
            set(gcf,'Position',[200,100,600,600])
%             set(gca,'FontSize',28,'FontName','Times New Roman')
            set(gca,'FontSize',40,'FontName','Arial')
            
            ax = gca;
            ax.Layer = 'top';
            
        end
    end
    toc
    pause(0.1)
end
% figure
% box on
%
% delta_th = 1e-3;
% p_n = patch(isosurface(XX, YY, ZZ, CC, -delta_th));
% p_p = patch(isosurface(XX, YY, ZZ, CC, delta_th));
%
% color_parula = parula(64);
%
% p_n.FaceColor = color_parula(1,:);
% p_p.FaceColor = color_parula(64,:);
% p_n.EdgeColor = 'none';
% p_p.EdgeColor = 'none';
%
% daspect([1 1 1])
% view([60 20]);
% axis tight
% camlight
% lighting gouraud
% material shiny
%
% xlim([-5 5])
% ylim([-5 5])
% zlim([-5 5])
%
% xlabel('x','FontSize',24,'Interpreter','latex')
% ylabel('y','FontSize',24,'Interpreter','latex')
% zlabel('z','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
% set(gcf,'Position',[200,100,600,600])
% set(gca,'FontSize',28,'FontName','Times New Roman')