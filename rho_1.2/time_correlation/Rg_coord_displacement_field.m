clear;
close all;

%% background
% triangle
theta = [0 1/3 2/3]*2*pi;
T = [cos(theta);sin(theta)];
pgon_T = polyshape(T(1,:),T(2,:));
T2 = [(T(:,1)+T(:,2))/2,T(:,2),[0;0]];
pgon_T2 = polyshape(T2(1,:),T2(2,:));

V_Sb = T(:,2);
V_Sc = T(:,1)+T(:,2)/2;

% V_Sb = V_Sb*sqrt(3/2);
% V_Sc = V_Sc*sqrt(3/2);

V_S6 = T(:,2)+T(:,1)/2;
V_S7 = T(:,1);

V_S6 = V_S6/norm(V_S6)*sqrt(3/2);
V_S7 = V_S7/norm(V_S7)*sqrt(3/2);

figure(1)
hold on
box on
pgon_BG = pgon_T2;
pT = plot(pgon_T,'FaceColor','#666666','LineWidth',2);
pBG = plot(pgon_BG,'FaceColor','w','LineWidth',2,'FaceAlpha',1);
daspect([1 1 1])
xlim([-1.5e-2 1.5e-2])
ylim([0 3e-2])
xticks(-0.04:0.01:0.04)
yticks(0:0.01:0.1)
xlabel('{\itS}_7','FontSize',28);
ylabel('{\itS}_6','FontSize',28);
set(gca,'FontSize',28,'FontName','Times New Roman')
set(gca,'LineWidth',2)
set(gca, 'Layer', 'top')
set([pT,pBG],'handlevisibility','off')
% set(gca,'Position',[0.14 0.17 0.54 0.77])
set(gcf,'Position',[2000,100,600,630])
ax = gca;
ax.XAxis.Exponent = -2;
ax.YAxis.Exponent = -2;

n_conf = 1;
n_trajectory = 5;
N = 6912;

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
ts = ts(ts<=1e3);

n_frame = length(ts);

t = ts;

iT = 0;

%% start loop
for T_i = fliplr([0.4]) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    %% load rg.mat file
    %     filename_rg = ['../rg_',num2str(T_i,'%.1f'),'.mat'];
    filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_rg)
    
    % assign color representing the temperature
    cp = parula(20);
    %     color = cp(round(T_i*10),:);
    
    for ic = 1 % loop over configuration
        disp(ic)
        
        %% time origin
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        [p_0,coord_0] = calculate_shape(rg_all);
        
        %%%%%%%%%%%%%%%% sort Rg into 3 type %%%%%%%%%%%%%%%%
        nstate = 3;
        state0 = zeros(1*N,1);
        [S_p,I_p] = sort(p_0(:,3));
        ni_s = floor(1*N/nstate);
        
        for is = 1
            state0(I_p(ni_s*(is-1)+1:ni_s*(is))) = is;
        end
        
        index_nonspere = find(state0==0);
        [S_S7,I_S7] = sort(coord_0(index_nonspere,2));
        for is = 2:3
            state0(index_nonspere(I_S7(ni_s*(is-2)+1:ni_s*(is-1)))) = is;
        end
        
        for i_s = 2
            index_shape = find(state0==i_s);
            index_rand = ceil(rand(n_trajectory,1)*length(index_shape));
            index_trajectory = index_shape(index_rand);
            
            for i_traj = 1:n_trajectory
                plot(coord_0(index_trajectory(i_traj),2),...
                    coord_0(index_trajectory(i_traj),1),'ko')
            end
            %% loop over time
            for it = 2:n_frame
                rg_all = rg_all_t_c(:,:,:,it,ic);
                [p_t,coord_t] = calculate_shape(rg_all);
                
                rg_all = rg_all_t_c(:,:,:,it-1,ic);
                [p_tm,coord_tm] = calculate_shape(rg_all);
                
                
                %% figures
                for i_traj = 1:n_trajectory
                    plot([coord_t(index_trajectory(i_traj),2),coord_tm(index_trajectory(i_traj),2)],...
                        [coord_t(index_trajectory(i_traj),1),coord_tm(index_trajectory(i_traj),1)],'k-')
                end
            end
            
        end
        
    end
    
    
    
    toc
end% end of temperature loop


