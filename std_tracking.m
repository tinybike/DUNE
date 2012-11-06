% std_tracking.m

max_t = max(tFinal);
prot1_k(isnan(prot1_k)) = 0;
prot2_k(isnan(prot2_k)) = 0;
prot3_k(isnan(prot3_k)) = 0;
prot4_k(isnan(prot4_k)) = 0;
prot5_k(isnan(prot5_k)) = 0;
prot6_k(isnan(prot6_k)) = 0;
prot11_k(isnan(prot11_k)) = 0;
prot101_k(isnan(prot101_k)) = 0;
prot1001_k(isnan(prot1001_k)) = 0;
for i = 1:max_t
    median_prot1_k(i) = median(nonzeros(prot1_k(:,i)));
    median_prot6_k(i) = median(nonzeros(prot6_k(:,i)));
    median_prot11_k(i) = median(nonzeros(prot11_k(:,i)));
    median_prot101_k(i) = median(nonzeros(prot101_k(:,i)));
    median_prot1001_k(i) = median(nonzeros(prot1001_k(:,i)));
end

prot1_b(isnan(prot1_b)) = 0;
prot2_b(isnan(prot2_b)) = 0;
prot3_b(isnan(prot3_b)) = 0;
prot4_b(isnan(prot4_b)) = 0;
prot5_b(isnan(prot5_b)) = 0;
prot6_b(isnan(prot6_b)) = 0;
prot11_b(isnan(prot11_b)) = 0;
prot101_b(isnan(prot101_b)) = 0;
prot1001_b(isnan(prot1001_b)) = 0;
size_of_b = size(prot1_b);
for i = 1:size_of_b(2)
    median_prot1_b(i) = median(nonzeros(prot1_b(:,i)));
    median_prot6_b(i) = median(nonzeros(prot6_b(:,i)));
    median_prot11_b(i) = median(nonzeros(prot11_b(:,i)));
    median_prot101_b(i) = median(nonzeros(prot101_b(:,i)));
    median_prot1001_b(i) = median(nonzeros(prot1001_b(:,i)));
end
