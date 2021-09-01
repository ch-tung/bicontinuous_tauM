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
index_time = [21];

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
            Ct_center = Ct-mean(Ct);
%             Ct_center = Ct-C_th;
            
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
            
            dl = Lx/100; % pixel size
            n_edge = 1;
            [XX, YY, ZZ] = meshgrid((-n_edge*dl:dl:Lx+n_edge*dl)-Lx/2,(-n_edge*dl:dl:Ly+n_edge*dl)-Ly/2,(-n_edge*dl:dl:Lz+n_edge*dl)-Ly/2);
            F = scatteredInterpolant(Xs, Ys, Zs, ps, 'natural');
            CC = F(XX, YY, ZZ);
%             CC = griddata(Xs, Ys, Zs, ps, XX, YY, ZZ, 'natural');
            IC = CC>0;
            f_C = sum(IC,'all')/length(XX(:)); % fraction of contious phase
            
            %% isosurface
            [f, v] = isosurface(XX, YY, ZZ, CC, 0);
            TR = triangulation(f,v);
%             [GC, MC, Area] = curvatures(v(:,1),v(:,2),v(:,3),f);
%             
%             i_interior = sum(min(abs(v),[Lx/2,Lx/2,Lx/2])==[Lx/2,Lx/2,Lx/2],2)==0;
%             
%             A_i = Area(i_interior);
%             GC_i = GC(i_interior);
%             integral_A = sum(A_i);
%             mean_GC = sum(GC_i.*A_i')/integral_A;
            
            %% figure
            figure
            box on
            
            boxsize = 10/2;
            boxcenter = [0 0 0];
            
            daspect([1 1 1])
            view([30 20])
%             set(gca, 'Projection','perspective')
           
%             TR = triangulation(f,v);
            Norms = faceNormal(TR);
            BC = boxcenter;
            back_facing = sum(Norms.*bsxfun(@minus,BC,campos),2)<=0;
            
            p = patch('Faces',f,'Vertices',v);
           
            color_parula = parula(64);
            
            load('YlOrRd.mat')
            color = rgba(:,1:3);
            color = gray(256);
            
%             color_dark = [139,0,0]/256;
%             color_bright = [255,192,64]/256;
            
            color_dark = color(64,:);
            color_bright = color(256,:);
            
%             color_dark = [0 0.16 0.64];
%             color_bright = [1 1 1];

%             color_dark = color_parula(16,:);
%             color_bright = color_parula(64,:);
            
            p.EdgeColor = 'none';
            p.FaceColor = 'flat';
            p.FaceVertexCData = repmat(color_bright,size(f,1),1).*back_facing +... 
                repmat(color_dark,size(f,1),1).*(1-back_facing);
            
            axis tight
            camlight('headlight')
            lighting gouraud    
            material dull 
            
            xlim([-boxsize boxsize]+boxcenter(1))
            ylim([-boxsize boxsize]+boxcenter(1))
            zlim([-boxsize boxsize]+boxcenter(1))
            
            xlabel('x','FontSize',24)
            ylabel('y','FontSize',24)
            zlabel('z','FontSize',24)
            
            set(gca,'LineWidth',2)
            set(gcf,'Position',[200,100,600,600])
            
            set(gca,'FontSize',40,'FontName','Arial')
            
            ax = gca;
            ax.Layer = 'top';
            
        end
    end
    toc
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