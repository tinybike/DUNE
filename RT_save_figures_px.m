% RT_save_figures_px.m

figure(1)
clf

offset = 1;

subplot(3,1,1)
load new_hsa_2
for i = 1:optimax
    semilogy(close_x,offset^2*clsy_mat(:,i)/sum(clsy_mat(:,i)),'Color',[0.75 1 0.75]);
    hold on
end
semilogy(data_close_x,offset^2*data_close_y,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10)
semilogy(close_x,offset^2*close_y,'g-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Closeness centrality','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
for i = 1:optimax
    semilogy(close_x,offset*clsy_mat(:,i)/sum(clsy_mat(:,i)),'Color',[0.75 0.75 1]);
    hold on
end
semilogy(data_close_x,offset*data_close_y,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
semilogy(close_x,offset*close_y,'b-','LineWidth',2)
% set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('$p(\ell)$','interpreter','latex','fontsize',22)

subplot(3,1,3)
load new_dme_1
for i = 1:optimax
    semilogy(close_x,clsy_mat(:,i)/sum(clsy_mat(:,i)),'Color',[1 0.75 0.75]);
    hold on
end
semilogy(data_close_x,data_close_y,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10)
semilogy(close_x,close_y,'r-','LineWidth',2)
set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.12 0.8 0.23])
xlabel('$\ell$','interpreter','latex','fontsize',22)

print('-depsc2','RT_px');