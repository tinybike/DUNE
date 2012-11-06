% RT_save_figures_frag.m

figure
clf

subplot(3,1,1)
load new_hsa_2
for i = 1:optimax
    plot(L1_del,L1_frag(i,:),'Color',[0.75 1 0.75]);
    hold on
    plot(L1_del,L1_pref_frag(i,:),'--','Color',[0.75 1 0.75]);
end
plot(L1_del,mean_L1_frag,'g-','LineWidth',2)
plot(data_L1_del,data_L1_frag,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10);
plot(L1_del,mean_L1_pref_frag,'g--');
plot(data_L1_del,data_L1_pref_frag,'gs','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10);
axis([-0.01 1.01 -0.1 1.1])
set(gca,'xtick',[])
set(gca,'Position',[0.12 0.64 0.8 0.25])
title('Error tolerance','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
for i = 1:optimax
    plot(L1_del,L1_frag(i,:),'Color',[0.75 0.75 1]);
    hold on
    plot(L1_del,L1_pref_frag(i,:),'--','Color',[0.75 0.75 1]);
end
plot(L1_del,mean_L1_frag,'b-','LineWidth',2)
plot(data_L1_del,data_L1_frag,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10);
plot(L1_del,mean_L1_pref_frag,'b--');
plot(data_L1_del,data_L1_pref_frag,'bs','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10);
axis([-0.01 1.01 -0.1 1.1])
ylabel('fragmentation','interpreter','latex','fontsize',22)
set(gca,'xtick',[])
set(gca,'Position',[0.12 0.38 0.8 0.25])

subplot(3,1,3)
load new_dme_1
for i = 1:optimax
    plot(L1_del,L1_frag(i,:),'Color',[1 0.75 0.75]);
    hold on
    plot(L1_del,L1_pref_frag(i,:),'--','Color',[1 0.75 0.75]);
end
plot(L1_del,mean_L1_frag,'r-','LineWidth',2)
plot(data_L1_del,data_L1_frag,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10);
plot(L1_del,mean_L1_pref_frag,'r--');
plot(data_L1_del,data_L1_pref_frag,'rs','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10);
axis([-0.01 1.01 -0.1 1.1])
set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.12 0.8 0.25])

xlabel('fraction of nodes deleted','interpreter','latex','fontsize',22)
print('-depsc2','RT_frag');
