% RT_save_figures_plambda.m

figure(1)
clf

offset = 1;

subplot(3,1,1)
load new_hsa_2
for i = 1:optimax
    semilogy(ex,offset^2*ey_mat(:,i)/sum(ey_mat(:,i)),'Color',[0.75 1 0.75]);
    hold on
end
loglog(data_ex,offset^2*data_ey,'go','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0],'MarkerSize',10)
loglog(ex,offset^2*ey,'g-','LineWidth',2)
set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.68 0.8 0.23])
title('Eigenvalue spectrum','interpreter','latex','fontsize',22)

subplot(3,1,2)
load new_sce_11
for i = 1:optimax
    semilogy(ex,offset*ey_mat(:,i)/sum(ey_mat(:,i)),'Color',[0.75 0.75 1]);
    hold on
end
loglog(data_ex,offset*data_ey,'bo','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
loglog(ex,offset*ey,'b-','LineWidth',2)
set(gca,'xtick',[])
% set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.40 0.8 0.23])
ylabel('$p(\lambda)$','interpreter','latex','fontsize',22)

subplot(3,1,3)
load new_dme_1
for i = 1:optimax
    semilogy(ex,ey_mat(:,i)/sum(ey_mat(:,i)),'Color',[1 0.75 0.75]);
    hold on
end
loglog(data_ex,data_ey,'ro','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',10)
loglog(ex,ey,'r-','LineWidth',2)
set(gca,'xtickMode', 'auto')
set(gca,'Position',[0.12 0.12 0.8 0.23])
xlabel('$\lambda$','interpreter','latex','fontsize',22)

print('-depsc2','RT_plambda');