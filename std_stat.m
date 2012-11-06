% std_stat.m

y = median(y_mat');
norm_y = sum(y);
y = y/norm_y;

cy = median(cy_mat');
norm_cy = sum(cy);
cy = cy/norm_cy;

by = median(by_mat');
norm_by = sum(by);
by = by/norm_by;

close_y = median(clsy_mat');
norm_clsy = sum(close_y);
close_y = close_y/norm_clsy;
close_x = clsx;

ey = median(ey_mat');
norm_ey = sum(ey);
ey = ey/norm_ey;

nvy = median(nvy_mat');
norm_nvy = sum(nvy);
nvy = nvy/norm_nvy;

knny = median(knny_mat');
norm_knny = sum(knny);
knny = knny/norm_knny;

for unindex = 1:maxUniqueK
    temp = mean_C_mat(unindex,:);
    mean_C(unindex) = median(nonzeros(temp));
end

for unindex = 1:maxUniqueK
    temp = mean_b_mat(unindex,:);
    mean_b(unindex) = median(nonzeros(temp));
end

for unindex = 1:maxUniqueK
    temp = mean_knn_mat(unindex,:);
    mean_knn(unindex) = median(nonzeros(temp));
end

for unindex = 1:maxUniqueK
    temp = mean_cls_mat(unindex,:);
    mean_cls(unindex) = median(nonzeros(temp));
end

uniq_GC_K = unique_k;

mean_L1_frag = median(L1_frag);
mean_L1_pref_frag = median(L1_pref_frag);
