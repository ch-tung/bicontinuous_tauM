clear
close all

n_particle = 6912;
N = n_particle;       % number of particles

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
ic = 1;

%% start loop
for T_i = fliplr([0.4]) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
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
    for it = 2:9:65
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
%         Ct_center = Ct-C_th;
        
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
        
        dl = 0.25; % pixel size
        [XX, YY, ZZ] = meshgrid((0:dl:Lx)-Lx/2,(0:dl:Ly)-Ly/2,0);
        CC = griddata(Xs, Ys, Zs, ps, XX, YY, ZZ, 'natural');
        
        figure(1)
        subplot(2,4,itime)
        [M,c] = contourf(XX, YY, CC, [-100 0 100]);
        
        %         image([0 l(1,2)-l(1,1)],[0 l(2,2)-l(2,1)],CC_rgb8)
        set(gca,'YDir','normal')
%         colormap(gray)
%         caxis([-1e-10 1e-10])
        
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
        
        set(gca,'LineWidth',2)
        set(gcf,'Position',[200,100,1200,600])
        set(gca,'FontSize',20,'FontName','Times New Roman')
        
    end
    
end
