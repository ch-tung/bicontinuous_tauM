clear
close all

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

t = ts/1000;

iT = 0;
index_t = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];

%% space
dk = 2*pi/100;
kk = 0:dk:2*pi;

%% figure
[XX,YY] = meshgrid(kk,t(index_t));

nu = 1;

C_kt = exp(-nu.*XX.^2.*YY);
pcolor(XX,YY,C_kt)
shading interp
set(gca, 'YScale', 'log')

caxis([0 1])
xlim([0 3])
yticks(10.^([-3,0,3]))

xlabel('k','FontSize',24)
ylabel('t','FontSize',24)
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')


