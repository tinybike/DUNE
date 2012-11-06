% RT_save_figures_Q.m

figure(1)
clf

subplot(3,1,1)
load new_hsa_alltimesteps_1
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
    plot(T(i,50:50:round(nnz(numMeasures(i,:)))),Q(i,1:nnz(numMeasures100(i,:))),'Color',[0.75 1 0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn(Q,tau,25,tFinal);
plot(xcoord,median_ycoord,'g-','LineWidth',2)
plot(xcoord,data_Q*ones(length(xcoord),1),'g--','LineWidth',2)
% ylim([0 1])
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Modularity','interpreter','latex','fontsize',22)

subplot(3,1,2)
clearvars
load new_sce_alltimesteps_1
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
%     plot(T(i,100:100:round(nnz(numMeasures(i,:)))),Q(i,1:nnz(numMeasures100(i,:))),'Color',[0.75 0.75 1]);
    plot(T(i,1:nnz(numMeasures(i,:))),Q(i,1:nnz(numMeasures(i,:))),'Color',[0.75,0.75,1]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(Q,tau,10,tFinal);
plot(xcoord,median_ycoord,'b-','LineWidth',2)
plot(xcoord,data_Q*ones(length(xcoord),1),'b--','LineWidth',2)
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('$Q$','interpreter','latex','fontsize',22)

subplot(3,1,3)
clearvars
load new_dme_alltimesteps_noeigs_2
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
	plot(T(i,1:nnz(numMeasures(i,:))),Q(i,1:nnz(numMeasures(i,:))),'Color',[1,0.75,0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(Q,tau,25,tFinal);
plot(xcoord,median_ycoord,'r-','LineWidth',2)
plot(xcoord,data_Q*ones(length(xcoord),1),'r--','LineWidth',2)
set(gca,'Position',[0.12 0.12 0.8 0.23])

xlabel('$t$ (Myr)','interpreter','latex','fontsize',22)
% ylim([0 1])
print('-depsc2','RT_Q');