% RT_save_figures_pb.m

figure(1)
clf

offset = 1;

subplot(3,1,1)
load new_hsa_2
for i = 1:optimax
    loglog(bx,offset^2*by_mat(:,i)/sum(by_mat(:,i)),'Color',[0.75 1 0.75]);
    hold on
end
loglog(data_bx,offset^2*data_by,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10)
loglog(bx,offset^2*by,'g-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Betweenness centrality','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
for i = 1:optimax
    loglog(bx,offset*by_mat(:,i)/sum(by_mat(:,i)),'Color',[0.75 0.75 1]);
    hold on
end
loglog(data_bx,offset*data_by,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
loglog(bx,offset*by,'b-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('$p(b)$','interpreter','latex','fontsize',22)

subplot(3,1,3)
load new_dme_1
for i = 1:optimax
    loglog(bx,by_mat(:,i)/sum(by_mat(:,i)),'Color',[1 0.75 0.75]);
    hold on
end
loglog(data_bx,data_by,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10)
loglog(bx,by,'r-','LineWidth',2)
set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.12 0.8 0.23])
xlabel('$b$','interpreter','latex','fontsize',22)

print('-depsc2','RT_pb');