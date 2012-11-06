% RT_save_figures_pk.m

figure(1)
clf

offset = 1;

subplot(3,1,1)
load new_hsa_2
for i = 1:optimax
    loglog(x,offset^2*y_mat(:,i)/sum(y_mat(:,i)),'Color',[0.75 1 0.75]);
    hold on
end
loglog(data_x,offset^2*data_y,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10)
loglog(x,offset^2*y,'g-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Degree centrality','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
for i = 1:optimax
    loglog(x,offset*y_mat(:,i)/sum(y_mat(:,i)),'Color',[0.75 0.75 1]);
    hold on
end
loglog(data_x,offset*data_y,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
loglog(x,offset*y,'b-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('$p(k)$','interpreter','latex','fontsize',22)

subplot(3,1,3)
load new_dme_1
for i = 1:optimax
    loglog(x,y_mat(:,i)/sum(y_mat(:,i)),'Color',[1 0.75 0.75]);
    hold on
end
loglog(data_x,data_y,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10)
loglog(x,y,'r-','LineWidth',2)
set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.12 0.8 0.23])

% Labels etc.
xlabel('$k$','interpreter','latex','fontsize',22)
print('-depsc2','RT_pk');