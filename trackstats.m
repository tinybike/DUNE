% trackstats.m

prot1_k(opt_m,t) = sum(A(:,1));
if N > 1
    prot2_k(opt_m,t) = sum(A(:,2));
end
if N > 2
    prot3_k(opt_m,t) = sum(A(:,3));
end
if N > 3
    prot4_k(opt_m,t) = sum(A(:,4));
end
if N > 4
    prot5_k(opt_m,t) = sum(A(:,5));
end
if N > 5
    prot6_k(opt_m,t) = sum(A(:,6));
end
if N > 10
    prot11_k(opt_m,t) = sum(A(:,11));
end
if N > 100
    prot101_k(opt_m,t) = sum(A(:,101));
end
if N > 1000
    prot1001_k(opt_m,t) = sum(A(:,1001));
end

b = 2/((N-1)*(N-2))*betweenness_centrality(A,'unweighted',1);
prot1_b(opt_m,t) = b(1);
if N > 1
    prot2_b(opt_m,t) = b(2);
end
if N > 2
    prot3_b(opt_m,t) = b(3);
end
if N > 3
    prot4_b(opt_m,t) = b(4);
end
if N > 4
    prot5_b(opt_m,t) = b(5);
end
if N > 5
    prot6_b(opt_m,t) = b(6);
end
if N > 10
    prot11_b(opt_m,t) = b(11);
end
if N > 100
    prot101_b(opt_m,t) = b(101);
end
if N > 1000
    prot1001_b(opt_m,t) = b(1001);
end