% RT_save_figures_loss.m

figure(1)
clf

subplot(3,1,1)
load new_hsa_disconnect_alltimesteps_1
for opt_m = 1:50
    orphanfrac(opt_m,:) = orphan(opt_m,:)./(nonorphan(opt_m,:) + orphan(opt_m,:));
end
gcFrac = orphanfrac;
gcFrac_le10 = orphanfrac;
gcFrac_nz = orphanfrac;
HBLC = orphanfrac;
Q = orphanfrac;
gcc = orphanfrac;
E2 = orphanfrac;
diameter = orphanfrac;
avgCls = orphanfrac;
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
    plot(T(i,1:nnz(numMeasures(i,:))),orphanfrac(i,1:nnz(numMeasures(i,:))),'Color',[0.75,1,0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(orphanfrac,tau,25,tFinal);
plot(xcoord,median_ycoord,'g-','LineWidth',2)
% plot(xcoord,data_Q*ones(length(xcoord),1),'g--','LineWidth',2)
%ylim([0 1])
set(gca,'Position',[0.12 0.68 0.8 0.23])
% title('Fraction of orphan proteins','interpreter','latex','fontsize',22)

subplot(3,1,2)
clearvars
load new_sce_disconnect_alltimesteps_1
for opt_m = 1:50
    orphanfrac(opt_m,:) = orphan(opt_m,:)./(nonorphan(opt_m,:) + orphan(opt_m,:));
end
gcFrac = orphanfrac;
gcFrac_le10 = orphanfrac;
gcFrac_nz = orphanfrac;
HBLC = orphanfrac;
Q = orphanfrac;
gcc = orphanfrac;
E2 = orphanfrac;
diameter = orphanfrac;
avgCls = orphanfrac;
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
%     plot(T(i,100:100:round(nnz(numMeasures(i,:)))),Q(i,1:nnz(numMeasures100(i,:))),'Color',[0.75 0.75 1]);
    plot(T(i,1:nnz(numMeasures(i,:))),orphanfrac(i,1:nnz(numMeasures(i,:))),'Color',[0.75,0.75,1]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(orphanfrac,tau,10,tFinal);
plot(xcoord,median_ycoord,'b-','LineWidth',2)
% plot(xcoord,data_Q*ones(length(xcoord),1),'b--','LineWidth',2)
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('fraction of orphan proteins','interpreter','latex','fontsize',22)

subplot(3,1,3)
clearvars
load new_dme_disconnect_alltimesteps_1
for opt_m = 1:50
    orphanfrac(opt_m,:) = orphan(opt_m,:)./(nonorphan(opt_m,:) + orphan(opt_m,:));
end
gcFrac = orphanfrac;
gcFrac_le10 = orphanfrac;
gcFrac_nz = orphanfrac;
HBLC = orphanfrac;
Q = orphanfrac;
gcc = orphanfrac;
E2 = orphanfrac;
diameter = orphanfrac;
avgCls = orphanfrac;
std_dyn;
T = cumsum(tau')';
for i = 1:optimax
	plot(T(i,1:nnz(numMeasures(i,:))),orphanfrac(i,1:nnz(numMeasures(i,:))),'Color',[1,0.75,0.75]);
    hold on
end
[xcoord,ycoord,median_ycoord] = ppi_dyn_ETS(orphanfrac,tau,25,tFinal);
plot(xcoord,median_ycoord,'r-','LineWidth',2)
% plot(xcoord,data_Q*ones(length(xcoord),1),'r--','LineWidth',2)
set(gca,'Position',[0.12 0.12 0.8 0.23])

xlabel('$t$ (Myr)','interpreter','latex','fontsize',22)
%ylim([0 1])
print('-depsc2','RT_loss');