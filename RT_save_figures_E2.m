% RT_save_figures_E2.m

figure(1)
clf

subplot(3,1,1)
load new_hsa_alltimesteps_1
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
    plot(T(i,50:50:round(nnz(numMeasures(i,:)))),E2(i,1:nnz(numMeasures100(i,:))),'Color',[0.75 1 0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn(E2,tau,25,tFinal);
plot(xcoord,median_ycoord,'g-','LineWidth',2)
[dmax_pos,dmax_pos] = max(data_E);
data_E(dmax_pos) = [];
[dmax_pos,dmax_pos] = max(data_E);
data_E2 = data_E(dmax_pos);
plot(xcoord,data_E2*ones(length(xcoord),1),'g--','LineWidth',2)
% ylim([0 1])
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Second-largest eigenvalue','interpreter','latex','fontsize',22)

subplot(3,1,2)
clearvars
load new_sce_alltimesteps_1
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
%     plot(T(i,100:100:round(nnz(numMeasures(i,:)))),E2(i,1:nnz(numMeasures100(i,:))),'Color',[0.75 0.75 1]);
    plot(T(i,1:nnz(numMeasures(i,:))),E2(i,1:nnz(numMeasures(i,:))),'Color',[0.75,0.75,1]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(E2,tau,10,tFinal);
plot(xcoord,median_ycoord,'b-','LineWidth',2)
[dmax_pos,dmax_pos] = max(data_E);
data_E(dmax_pos) = [];
[dmax_pos,dmax_pos] = max(data_E);
data_E2 = data_E(dmax_pos);
plot(xcoord,data_E2*ones(length(xcoord),1),'b--','LineWidth',2)
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylim([-1 1])
ylabel('$\lambda_2$','interpreter','latex','fontsize',22)

subplot(3,1,3)
clearvars
load new_dme_alltimesteps_noeigs_2
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
	plot(T(i,1:nnz(numMeasures(i,:))),E2(i,1:nnz(numMeasures(i,:))),'Color',[1,0.75,0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(E2,tau,25,tFinal);
plot(xcoord,median_ycoord,'r-','LineWidth',2)
[dmax_pos,dmax_pos] = max(data_E);
data_E(dmax_pos) = [];
[dmax_pos,dmax_pos] = max(data_E);
data_E2 = data_E(dmax_pos);
plot(xcoord,data_E2*ones(length(xcoord),1),'r--','LineWidth',2)
ylim([-1 1])
set(gca,'Position',[0.12 0.12 0.8 0.23])

xlabel('$t$ (Myr)','interpreter','latex','fontsize',22)
% ylim([0 1])
print('-depsc2','RT_E2');