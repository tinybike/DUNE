% ppi_import.m

if ~exist('data_C','var')
    fprintf('Importing data...');
%     data_B = sparse(importdata('hitpredict/hsa_ss_hiconf_A.csv',','));
%     data_B = sparse(importdata('hitpredict/sce_ss_hiconf_A.csv',','));
    data_B = sparse(importdata('hitpredict/dme_ss_hiconf_A.csv',','));
%     data_B = sparse(importdata('hitpredict/eco_ss_hiconf_A.csv',','));
%     data_B = sparse(importdata('hitpredict/eco_predicted_A.csv',','));
%     data_B = sparse(importdata('bacteriome/eco_functional_A.csv',','));
%     data_B = sparse(importdata('hitpredict/dme_predicted_A.csv',','));
%     data_B = sparse(importdata('hitpredict/rat_ss_hiconf_A.csv',','));
%     data_B = sparse(importdata('hitpredict/mus_ss_hiconf_A.csv',','));
%     data_B = sparse(importdata('hitpredict/spo_ss_hiconf_A.csv',','));
    data_disconnect = find(sum(data_B) == 0);
    for data_disc_i = 1:length(data_disconnect)
        data_B(:,1:length(data_B) == data_disconnect(data_disc_i) - data_disc_i + 1) = [];
        data_B(1:length(data_B) == data_disconnect(data_disc_i) - data_disc_i + 1,:) = [];
    end
    cutoff = round(1*length(data_B));
    data_A = data_B(1:cutoff,1:cutoff);
    fprintf('calculating statistics...');
    % Full network degree
    data_results = full(sum(data_A)');
    % L_1 degree only
    [data_ci data_sizes] = components(data_A);
    [data_max_comp_size,data_max_comp] = max(data_sizes);
    data_compA = data_A;
%     data_results = full(sum(data_compA)');
    mean_data_results = mean(data_results);
    data_compA(:,data_ci ~= data_max_comp) = [];
    data_compA(data_ci ~= data_max_comp,:) = [];
    data_dist = all_shortest_paths(data_compA);
    data_cls = sum(data_dist)/(data_max_comp_size - 1);
    data_cls = 1./data_cls;
    data_b = 2/((data_max_comp_size-1)*(data_max_comp_size-2))*betweenness_centrality(data_compA,'unweighted',1);
    %data_b = betweenness_centrality(data_compA,'unweighted',1);
    %data_b = data_b/sum(data_b);
    data_GC_results = sum(data_compA);
    [data_EV,data_E] = eig(full(data_compA)./repmat(sum(data_compA),data_max_comp_size,1));
    data_E = real(sum(data_E));
    [data_max_pos,data_max_pos] = max(data_E);
    data_NV = data_EV(:,data_max_pos);
    clear data_EV
    data_NV = data_NV/sum(data_NV);
    data_C = clustering_coefficients(data_compA,'unweighted',1);
    data_GC_C = data_C;
    data_knn = zeros(data_max_comp_size,1);
    for i = 1:data_max_comp_size
        data_nzA = find(data_compA(:,i));
        data_knn(i) = median(data_GC_results(data_nzA));
    end
    %data_knn = data_knn/(data_max_comp_size-1);
    data_N = nnz(sum(data_compA));
    %data_N = nnz(data_A);
    data_bins = max(data_results);
    data_edges = full(sum(sum(data_compA)) + trace(data_compA))/2;
    data_nd = 1./data_cls/max(1./data_cls);
    % Fragmentation
    data_L1 = data_compA;
    data_L1_pref = data_compA;
    data_L1_frag = zeros(1,frags);
    data_L1_pref_frag = zeros(1,frags);
    data_frac_del = (1:data_max_comp_size)/data_max_comp_size;
    data_L1_del = linspace(0,1,frags);
    for i = 1:frags
        data_num_del(i) = nnz(data_frac_del < data_L1_del(i));
        if data_num_del(i)
            data_frac_del(1:data_num_del(i)) = [];
        end
        for j = 1:data_num_del(i)
            [ndel,ndel] = max(rand(1,length(data_L1)));
            [pdel,pdel] = max(sum(data_L1_pref));
            data_L1(ndel,:) = [];
            data_L1(:,ndel) = [];
            data_L1_pref(pdel,:) = [];
            data_L1_pref(:,pdel) = [];
        end
        [data_L1_ci data_L1_sizes] = components(data_L1);
        data_L1_frag(i) = max(data_L1_sizes)/data_max_comp_size;
        [data_L1_pref_ci data_L1_pref_sizes] = components(data_L1_pref);
        data_L1_pref_frag(i) = max(data_L1_pref_sizes)/data_max_comp_size;
    end
    % frag10
    data_L1 = data_compA;
    data_frac_del = (1:data_max_comp_size)/data_max_comp_size;
    data_del10 = nnz(data_frac_del < fragfrac);
    for j = 1:data_del10
        [ndel,ndel] = max(rand(1,length(data_L1)));
        data_L1(ndel,:) = [];
        data_L1(:,ndel) = [];
    end
    [data_L1_ci data_L1_sizes] = components(data_L1);
    data_frag10 = max(data_L1_sizes)/data_max_comp_size;
    % HBLC
    data_k = full(sum(data_compA))';
    data_HBLC = nnz((data_k < mean(data_k)).*(data_b > mean(data_b)))/data_max_comp_size;
    [data_Ci data_Q] = modularity_louvain_und(data_compA);
    data_Eglob = efficiency(data_compA);
    data_Eloc = efficiency(data_compA,1);
    % calculate mean overlays for scatter plots
    data_GC_uniq_K = unique(data_GC_results');
    count_GC_X = length(data_GC_uniq_K);
    data_mean_C = zeros(count_GC_X,1);
    data_mean_b = zeros(count_GC_X,1);
    data_mean_knn = zeros(count_GC_X,1);
    data_mean_cls = zeros(count_GC_X,1);
    for i = 1:count_GC_X
        find_values = find(data_GC_results == i);
        if sum(find_values)
            data_mean_C(i) = median(data_C(find_values));
            data_mean_b(i) = median(data_b(find_values));
            data_mean_knn(i) = median(data_knn(find_values));
            data_mean_cls(i) = median(data_cls(find_values));
        end
%         find_values = (data_GC_results' == data_GC_uniq_K(i)).*data_C;
%         find_values(find_values == 0) = [];
%         data_mean_C(i) = median(find_values);
%         find_values = (data_GC_results' == data_GC_uniq_K(i)).*data_b;
%         find_values(find_values == 0) = [];
%         data_mean_b(i) = median(find_values);
%         find_values = (data_GC_results' == data_GC_uniq_K(i)).*data_knn;
%         find_values(find_values == 0) = [];
%         data_mean_knn(i) = median(find_values);
%         find_values = (data_GC_results == data_GC_uniq_K(i)).*data_cls;
%         find_values(find_values == 0) = [];
%         data_mean_cls(i) = median(find_values);
    end
    % Histogram experimental data
    avgBinsResults = linspace(1,100,numBinsResults);
    avgBins = linspace(0,1,numBins);
    avgBinsB = linspace(0,max(data_b)*1.5,numBins);
    avgBinsCls = linspace(0,max(data_cls)*2,numBins);
    avgBinsEigs = linspace(-1,1,numBinsResults);
    avgBinsNV = linspace(0,max(data_NV)*1.5,numBins);
    [data_y,data_x] = hist(data_results,avgBinsResults);
    data_y = data_y/sum(data_y);
    [data_ey,data_ex] = hist(data_E,avgBinsEigs);
    data_ey = data_ey/sum(data_ey);
    [data_cy,data_cx] = hist(data_C,avgBins);
    data_cy = data_cy/sum(data_cy);
    [data_by,data_bx] = hist(data_b,avgBinsB);
    data_by = data_by/sum(data_by);
    [data_close_y,data_close_x] = hist(data_cls,avgBinsCls);
    data_close_y = data_close_y/sum(data_close_y);
    [data_nvy,data_nvx] = hist(data_NV,avgBinsNV);
    data_nvy = data_nvy/sum(data_nvy);
    [data_ndy,data_ndx] = hist(data_nd,avgBins);
    data_ndy = data_ndy/sum(data_ndy);
    [data_knny,data_knnx] = hist(data_knn,avgBinsResults);
    data_knny = data_knny/sum(data_knny);
    fprintf('done.\n');
end
