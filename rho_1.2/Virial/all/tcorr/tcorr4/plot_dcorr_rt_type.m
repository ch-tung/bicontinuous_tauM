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

index_t = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];

sim_r_c = sim_r_10_c;
%     s_fluc = sqrt(sim_s_c(:,12));
s_fluc = 0*sqrt(1e-2);
sim_r_normalized = (sim_r_c-sim_a_c)./(sim_s_c-sim_a_c + s_fluc^2);

rr = 0:0.25:10;
sim_r_normalized = [sim_s_c; sim_r_normalized];


% figures
%% C(r,t)
figure
box on
log_time = log(ts(2:end));
[X,Y] = meshgrid(rr([1 5:41]),ts(index_t));
%     sim_r_normalized(sim_r_normalized<1e-8) = 1e-8;
pcolor(X,Y/1e3,sim_r_normalized([1 5:41],(index_t))')
%     pcolor(X,Y,(sim_r_c(:,(2:end))'))
hold on
%     contour(X,Y,(sim_r_normalized(:,(2:end))'),'w')
contour(X,Y/1e3,sim_r_normalized([1 5:41],(index_t))',0.1:0.1:1,'w')

%     for n = 0:7
%         v = 10^(-n);
%         plot(0:0.1:10,(0:0.1:10)/v,'--w')
%     end

shading interp
set(gca, 'YScale', 'log')
%     set(gca,'ColorScale','log')
xlim([1 5])
ylim([1e0 1e7]/1e3)
xticks([1:10])
yticks(10.^([-3,0,3]))
caxis([-0.5 0.5])

ax = gca;
ax.Layer = 'top';

%     colorbar

%     colormap(jet)

xlabel('$r$','FontSize',24,'Interpreter','latex')
ylabel('$t$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Times New Roman')
