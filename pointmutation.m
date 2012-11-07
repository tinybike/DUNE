function A = pointmutation(A,linker,c,max_neighbors)

N = length(A);
targetlink = rand(1,N);
[target,target] = max(targetlink);
neighbors = cell(max_neighbors);
subtarget = zeros(max_neighbors,1);
ind_target = zeros(max_neighbors,1);
if sum(A(target,:)) > 0
    neighbors{1} = find(A(target,:));
    [subtarget(1),subtarget(1)] = max(rand(length(neighbors{1}),1));
    if rand(1) < c^2
        ind_target(1) = neighbors{1}(subtarget(1));
    end
    for nni = 2:max_neighbors
        if sum(A(neighbors{nni-1}(subtarget(nni-1)),:)) > 0
            neighbors{nni} = find(A(neighbors{nni-1}(subtarget(nni-1)),:));
            [subtarget(nni),subtarget(nni)] = max(rand(length(neighbors(nni)),1));
            if rand(1) < c^(nni+1)
                ind_target(nni) = neighbors{nni}(subtarget(nni));
            end
        end
    end
end

A(linker,target) = 1;
A(target,linker) = 1;
for ind_i = 1:max_neighbors
    if ind_target(ind_i)
        A(linker,ind_target(ind_i)) = 1;
        A(ind_target(ind_i),linker) = 1;
    end
end

