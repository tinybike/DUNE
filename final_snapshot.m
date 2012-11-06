% Fragmentation
[ci sizes] = components(A);
[max_comp_size,max_comp] = max(sizes);
compA = A;
compA(:,ci ~= max_comp) = [];
compA(ci ~= max_comp,:) = [];
L1 = compA;
L1_pref = compA;
frac_del = (1:max_comp_size)/max_comp_size;
L1_del = linspace(0,1,frags);
for i = 1:frags
    num_del(i) = nnz(frac_del < L1_del(i));
    if num_del(i)
        frac_del(1:num_del(i)) = [];
    end
    for j = 1:num_del(i)
        [ndel,ndel] = max(rand(1,length(L1)));
        [pdel,pdel] = max(sum(L1_pref));
        L1(ndel,:) = [];
        L1(:,ndel) = [];
        L1_pref(pdel,:) = [];
        L1_pref(:,pdel) = [];
    end
    [L1_ci L1_sizes] = components(L1);
    L1_frag(opt_m,i) = max(L1_sizes)/max_comp_size;
    [L1_pref_ci L1_pref_sizes] = components(L1_pref);
    L1_pref_frag(opt_m,i) = max(L1_pref_sizes)/max_comp_size;
end

% Remove disconnected nodes for stats
disconnect = find(sum(A) == 0);
for disc_i = 1:length(disconnect)
    A(:,1:N == disconnect(disc_i) - disc_i + 1) = [];
    A(1:N == disconnect(disc_i) - disc_i + 1,:) = [];
end
N = length(A);
results = full(sum(A)');
if sum(results) == 0
    continue;
end
[ci sizes] = components(A);
[max_comp_size,max_comp] = max(sizes);
compA = A;
compA(:,ci ~= max_comp) = [];
compA(ci ~= max_comp,:) = [];
%     compNormA = compA./repmat(sum(compA),max_comp_size,1);
Kinv = diag(sparse(1./sum(compA)));
compNormA = compA*Kinv;
b = 2/((max_comp_size-1)*(max_comp_size-2))*betweenness_centrality(compA,'unweighted',1);
C = clustering_coefficients(compA,'unweighted',1);
GC_results = full(sum(compA));
knn = zeros(max_comp_size,1);
for i = 1:max_comp_size
    nzA = find(compA(:,i));
    if isempty(nzA)
        nzA = 1;
    end
    knn(i) = mean(GC_results(nzA));
end
dist = all_shortest_paths(compA);
cls = sum(dist)/(max_comp_size-1);
cls = 1./cls;
%     [EV,E] = eig(full(compNormA));
E = real(eig(full(compNormA)));
%     if max_comp_size < 1000
%         E = real(eig(full(compNormA)));
%     else
%         E = real(eigs(compNormA,1000));
%     end
%     E = real(sum(E));
%     [max_pos,max_pos] = max(E);
%     NV = EV(:,max_pos);
%     clear EV
%     NV = NV/sum(NV);
% E = zeros(N,1);
NV = zeros(N,1);
%     Eloc = efficiency(compA,1);

nd = 1./cls/max(1./cls);

[y,x] = hist(results,avgBinsResults);
[cy,cx] = hist(C,avgBins);
[by,bx] = hist(b,avgBinsB);
[clsy,clsx] = hist(cls,avgBinsCls);
[ey,ex] = hist(E,avgBinsEigs);
[nvy,nvx] = hist(NV,avgBinsNV);
[knny,knnx] = hist(knn,avgBinsResults);
% [ndy,ndx] = hist(nd,avgBins);
[Ci Q_final] = modularity_louvain_und(compA);
%     [Qy,Qx] = hist(Ci,avgBinsQ);
%     [Elocy,Elocx] = hist(Eloc,avgBinsEigs);

unique_k = 1:maxUniqueK;
median_C = zeros(maxUniqueK,1);
median_b = zeros(maxUniqueK,1);
median_knn = zeros(maxUniqueK,1);
median_cls = zeros(maxUniqueK,1);
for i = 1:maxUniqueK
    find_values = find(GC_results == i);
    if sum(find_values)
        median_C(i) = median(C(find_values));
        median_b(i) = median(b(find_values));
        median_knn(i) = median(knn(find_values));
        median_cls(i) = median(cls(find_values));
    end
    
%     find_values = (GC_results(:) == i).*C;
%     if sum(find_values)
%         %             mean_C(i) = sum(find_values)/nnz(find_values);
%         median_C(i) = median(nonzeros(find_values));
%     end
%     find_values = (GC_results(:) == i).*b;
%     if sum(find_values)
%         %             mean_b(i) = sum(find_values)/nnz(find_values);
%         median_b(i) = median(nonzeros(find_values));
%     end
%     find_values = (GC_results(:) == i).*knn;
%     if sum(find_values)
%         %             mean_knn(i) = sum(find_values)/nnz(find_values);
%         median_knn(i) = median(nonzeros(find_values));
%     end
%     find_values = (GC_results == i).*cls;
%     if sum(find_values)
%         %             mean_cls(i) = sum(find_values)/nnz(find_values);
%         median_cls(i) = median(nonzeros(find_values));
%     end
end

total_mean_C = total_mean_C + median_C;
mean_C_mat(:,opt_m) = median_C';
total_mean_b = total_mean_b + median_b;
mean_b_mat(:,opt_m) = median_b';
total_mean_knn = total_mean_knn + median_knn;
mean_knn_mat(:,opt_m) = median_knn';
total_mean_cls = total_mean_cls + median_cls;
mean_cls_mat(:,opt_m) = median_cls';

total_y = total_y + y;
y_mat(:,opt_m) = y';
total_cy = total_cy + cy;
cy_mat(:,opt_m) = cy';
total_by = total_by + by;
by_mat(:,opt_m) = by';
total_clsy = total_clsy + clsy;
clsy_mat(:,opt_m) = clsy';
total_ey = total_ey + ey;
ey_mat(:,opt_m) = ey';
total_nvy = total_nvy + nvy;
nvy_mat(:,opt_m) = nvy';
total_knny = total_knny + knny;
knny_mat(:,opt_m) = knny';
% total_ndy = total_ndy + ndy;
% ndy_mat(:,opt_m) = ndy';
%     total_Qy = total_Qy + Qy;
%     Qy_mat(:,opt_m) = Qy';
%     total_Elocy = total_Elocy + Elocy;
%     Elocy_mat(:,opt_m) = Elocy';

diameter(opt_m) = max(dist(:));
[Ci Q(opt_m)] = modularity_louvain_und(compA);
f1(opt_m) = max_comp_size/N;
avgk(opt_m) = mean(results);
gcc(opt_m) = mean(C);
Nlist(opt_m) = N;
Klist(opt_m) = (sum(sum(A)) + trace(A))/2;
if exist('t','var')
tlist(opt_m) = t;
end

sizesList{opt_m} = sizes;

A_cell{opt_m} = A;
k_cell{opt_m} = results;
b_cell{opt_m} = b;
C_cell{opt_m} = C;
x_cell{opt_m} = cls;
disconnect_list(opt_m) = length(disconnect);
