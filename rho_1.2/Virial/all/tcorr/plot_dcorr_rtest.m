clear;
close all;

filetype = 'vp';
corrtype = {'11','10'};
covtype = 1;

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
index_t = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
% index_t = 2:2:length(t);
% index_t = 1:length(t);

%% start loop
% for T_i = ([0.4:0.1:1.0 1.2 1.5 2.0 3.0 5.0]) % loop over temperature
for T_i = 5.0 % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load dcorr_py.mat file
    filename_dcorr = ['tcorr4_test1.0','.mat'];
    load(filename_dcorr)
    
    for i_corr = 1:2
        sim_r_c = eval(['sim_r_',corrtype{i_corr},'_c']);
        
        % assign color representing the temperature
        cp = parula(20);
        
        %     s_fluc = sqrt(sim_s_c(:,12));
        s_fluc = 0*sqrt(1e-2);
        sim_r_normalized = (sim_r_c-sim_a_c)./(sim_s_c + s_fluc^2);
        
        if covtype==1
            sim_r_normalized = (sim_r_c);
        end
        
        
        % figures
        %% C(r)
        figure(iT)
        box on
        hold on

        plot(t,sim_r_normalized(10,:))
        set(gca, 'XScale', 'log')
%         ylim([-0.05 0.1])
        
        toc
    end
end% end of temperature loop
% close all;