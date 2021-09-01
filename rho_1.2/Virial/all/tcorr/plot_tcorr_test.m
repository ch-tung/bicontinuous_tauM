clear
close all


load('tcorr_CG_vo_py_0.5.mat')

%%
figure
hold on

% plot(t,C_0./Var_s)
% plot(t,C_1./Var_s)
plot(t,C_0./C_0(1))
plot(t,C_1./C_1(1))

set(gca, 'XScale', 'log')