% RT_save_figures_nvk.m

figure(1)
clf

subplot(3,1,1)
load new_hsa_2
std_stat;
for i = 1:optimax
    loglog(uniq_GC_K,mean_knn_mat(:,i),'Color',[0.75 1 0.75]);
    hold on
end
loglog(data_GC_uniq_K,data_mean_knn,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10)
hold on
loglog(uniq_GC_K,mean_knn,'g-','LineWidth',2)
ylim([1 50])
xlim([1 120])
% set(gca,'xtick',[])
% set(gca,'Position',[0.12 0.64 0.8 0.25])
title('Assortativity','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
std_stat;
for i = 1:optimax
    loglog(uniq_GC_K,mean_knn_mat(:,i),'Color',[0.75 0.75 1]);
    hold on
end
loglog(data_GC_uniq_K,data_mean_knn,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
hold on
loglog(uniq_GC_K,mean_knn,'b-','LineWidth',2);
ylim([1 50])
xlim([1 120])
% set(gca,'xtick',[])
% set(gca,'Position',[0.12 0.38 0.8 0.25])
ylabel('$\widetilde{n}$','interpreter','latex','fontsize',22)

subplot(3,1,3)
load new_dme_1
std_stat;
for i = 1:optimax
    loglog(uniq_GC_K,mean_knn_mat(:,i),'Color',[1 0.75 0.75]);
    hold on
end
loglog(data_GC_uniq_K,data_mean_knn,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10)
hold on
loglog(uniq_GC_K,mean_knn,'r-','LineWidth',2);
ylim([1 50])
xlim([1 120])
set(gca,'xtickMode', 'auto')
% set(gca,'Position',[0.12 0.12 0.8 0.25])
xlabel('$k$','interpreter','latex','fontsize',22)
print('-depsc2','RT_nvk');