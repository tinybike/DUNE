% std_dyn.m

% Set up a binary matrix for how long each simulation lasted
T = cumsum(tau')';
tFinal = zeros(1,optimax);
[maxNodes,maxNodes] = size(totalNumNodes);
measures = 1:maxNodes;
if length(gcFrac) ~= maxNodes
    measures100 = 50:50:maxNodes;
else
    measures100 = measures;
end
numMeasures = zeros(optimax,maxNodes);
[nM100 nM100] = size(gcFrac);
% nM100 = length(measures100);
numMeasures100 = zeros(optimax,nM100);
for i = 1:optimax
    tFinal(i) = max(totalNumNodes(i,:));
    numMeasures(i,:) = measures < tFinal(i);
    numMeasures100(i,:) = measures100 < tFinal(i);
end
fiftyPercent = sum(sum(numMeasures) > optimax/2);
fiftyPercent100 = sum(sum(numMeasures100) > optimax/2);

for iD = 1:maxNodes
    median_totalNumNodes(iD) = median(nonzeros(totalNumNodes(:,iD)));
end

% Mean + std of statistics calculated every time step
% Total number of edges
% mean_totalNumEdges = sum(totalNumEdges)./sum(numMeasures);
for iD = 1:maxNodes
    mean_totalNumEdges(iD) = median(nonzeros(totalNumEdges(:,iD)));
end
% std_totalNumEdges = sqrt(sum((totalNumEdges - numMeasures.*repmat(mean_totalNumEdges,optimax,1)).^2)./(sum(numMeasures) - 1));
% std_totalNumEdges(isinf(std_totalNumEdges)) = 0;
% std_totalNumEdges(isnan(std_totalNumEdges)) = 0;
% std_totalNumEdges(isinf(mean_totalNumEdges)) = [];
mean_totalNumEdges(isinf(mean_totalNumEdges)) = [];

% Fraction of random edges
% mean_fracRandomEdges = sum(fracRandomEdges)./sum(numMeasures);
fracRandomEdges(isnan(fracRandomEdges)) = 0;
for iD = 1:maxNodes
    mean_fracRandomEdges(iD) = median(nonzeros(fracRandomEdges(:,iD)));
end
% std_fracRandomEdges = sqrt(sum((fracRandomEdges - numMeasures.*repmat(mean_fracRandomEdges,optimax,1)).^2)./(sum(numMeasures) - 1));
% std_fracRandomEdges(isinf(std_fracRandomEdges)) = 0;
% std_fracRandomEdges(isnan(std_fracRandomEdges)) = 0;
% std_fracRandomEdges(isinf(mean_fracRandomEdges)) = [];
mean_fracRandomEdges(isinf(mean_fracRandomEdges)) = [];
% Fraction of duplicate edges
% mean_fracDuplicateEdges = sum(fracDuplicateEdges)./sum(numMeasures);
fracDuplicateEdges(isnan(fracDuplicateEdges)) = 0;
for iD = 1:maxNodes
    mean_fracDuplicateEdges(iD) = median(nonzeros(fracDuplicateEdges(:,iD)));
end
% std_fracDuplicateEdges = sqrt(sum((fracDuplicateEdges - numMeasures.*repmat(mean_fracDuplicateEdges,optimax,1)).^2)./(sum(numMeasures) - 1));
% std_fracDuplicateEdges(isinf(std_fracDuplicateEdges)) = 0;
% std_fracDuplicateEdges(isnan(std_fracDuplicateEdges)) = 0;
% std_fracDuplicateEdges(isinf(mean_fracDuplicateEdges)) = [];
mean_fracDuplicateEdges(isinf(mean_fracDuplicateEdges)) = [];

% Mean + std of statistics calculated every 100 time steps
% Global clustering coefficient
% mean_gcc = sum(gcc)./sum(numMeasures100);
for iD = 1:nM100
    mean_gcc(iD) = median(nonzeros(gcc(:,iD)));
end
% std_gcc = sqrt(sum((gcc - numMeasures100.*repmat(mean_gcc,optimax,1)).^2)./(sum(numMeasures100) - 1));
% std_gcc(isinf(std_gcc)) = 0;
% std_gcc(isnan(std_gcc)) = 0;
% std_gcc(isinf(mean_gcc)) = [];
mean_gcc(isinf(mean_gcc)) = [];
% Diameter
% mean_diameter = sum(diameter)./sum(numMeasures100);
diameter(isnan(diameter)) = 0;
diameter(isinf(diameter)) = 0;
for iD = 1:nM100
    mean_diameter(iD) = median(nonzeros(diameter(:,iD)));
end
% std_diameter = sqrt(sum((diameter - numMeasures100.*repmat(mean_diameter,optimax,1)).^2)./(sum(numMeasures100) - 1));
% std_diameter(isinf(std_diameter)) = 0;
% std_diameter(isnan(std_diameter)) = 0;
% std_diameter(isinf(mean_diameter)) = [];
mean_diameter(isinf(mean_diameter)) = [];
% Average distance
mean_avgCls = sum(avgCls)./sum(numMeasures100);
% std_avgCls = sqrt(sum((avgCls - numMeasures100.*repmat(mean_avgCls,optimax,1)).^2)./(sum(numMeasures100) - 1));
% std_avgCls(isinf(std_avgCls)) = 0;
% std_avgCls(isnan(std_avgCls)) = 0;
% std_avgCls(isinf(mean_avgCls)) = [];
mean_avgCls(isinf(mean_avgCls)) = [];
% Fraction of nodes in giant component
% mean_gcFrac = sum(gcFrac)./sum(numMeasures100);
for iD = 1:nM100
    mean_gcFrac(iD) = median(nonzeros(gcFrac(:,iD)));
end
% std_gcFrac = sqrt(sum((gcFrac - numMeasures100.*repmat(mean_gcFrac,optimax,1)).^2)./(sum(numMeasures100) - 1));
% std_gcFrac(isinf(std_gcFrac)) = 0;
% std_gcFrac(isnan(std_gcFrac)) = 0;
% std_gcFrac(isinf(mean_gcFrac)) = [];
mean_gcFrac(isinf(mean_gcFrac)) = [];
% Second-largest eigenvalue
% mean_E2 = sum(E2)./sum(numMeasures100);
E2(isnan(E2)) = 0;
E2(isinf(E2)) = 0;
for iD = 1:nM100
    mean_E2(iD) = median(nonzeros(E2(:,iD)));
end
% std_E2 = sqrt(sum((E2 - numMeasures100.*repmat(mean_E2,optimax,1)).^2)./(sum(numMeasures100) - 1));
% std_E2(isinf(std_E2)) = 0;
% std_E2(isnan(std_E2)) = 0;
% std_E2(isinf(mean_E2)) = [];
mean_E2(isinf(mean_E2)) = [];
if exist('frag10','var')
    % 10% fragmentation tolerance
    mean_frag10 = sum(frag10)./sum(numMeasures100);
%     std_frag10 = sqrt(sum((frag10 - numMeasures100.*repmat(mean_frag10,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_frag10(isinf(std_frag10)) = 0;
%     std_frag10(isnan(std_frag10)) = 0;
%     std_frag10(isinf(mean_frag10)) = [];
    mean_frag10(isinf(mean_frag10)) = [];
end
if exist('Q','var')
    % Fraction of nodes in giant component ignoring microcomponents
    mean_gcFrac_le10 = sum(gcFrac_le10)./sum(numMeasures100);
%     std_gcFrac_le10 = sqrt(sum((gcFrac_le10 - numMeasures100.*repmat(mean_gcFrac_le10,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_gcFrac_le10(isinf(std_gcFrac_le10)) = 0;
%     std_gcFrac_le10(isnan(std_gcFrac_le10)) = 0;
%     std_gcFrac_le10(isinf(mean_gcFrac_le10)) = 1;
    mean_gcFrac_le10(isinf(mean_gcFrac_le10)) = 1;
    % Fraction of nodes in giant component ignoring orphan nodes
%     mean_gcFrac_nz = sum(gcFrac_nz)./sum(numMeasures100);
    gcFrac_nz(isnan(gcFrac_nz)) = 0;
    gcFrac_nz(isinf(gcFrac_nz)) = 0;
    for iD = 1:nM100
        mean_gcFrac_nz(iD) = median(nonzeros(gcFrac_nz(:,iD)));
    end
%     std_gcFrac_nz = sqrt(sum((gcFrac_nz - numMeasures100.*repmat(mean_gcFrac_nz,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_gcFrac_nz(isinf(std_gcFrac_nz)) = 0;
%     std_gcFrac_nz(isnan(std_gcFrac_nz)) = 0;
%     std_gcFrac_nz(isinf(mean_gcFrac_nz)) = [];
    mean_gcFrac_nz(isinf(mean_gcFrac_nz)) = [];
    % Fraction of HBLC (high betweenness low connectivity) nodes
    mean_HBLC = sum(HBLC)./sum(numMeasures100);
%     std_HBLC = sqrt(sum((HBLC - numMeasures100.*repmat(mean_HBLC,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_HBLC(isinf(std_HBLC)) = 0;
%     std_HBLC(isnan(std_HBLC)) = 0;
%     std_HBLC(isinf(mean_HBLC)) = [];
    mean_HBLC(isinf(mean_HBLC)) = [];
    % Modularity
%     mean_Q = sum(Q)./sum(numMeasures100);
    Q(isnan(Q)) = 0;
    Q(isinf(Q)) = 0;
    for iQ = 1:nM100
        mean_Q(iQ) = median(nonzeros(Q(:,iQ)));
    end
%     std_Q = sqrt(sum((Q - numMeasures100.*repmat(mean_Q,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_Q(isinf(std_Q)) = 0;
%     std_Q(isnan(std_Q)) = 0;
%     std_Q(isinf(mean_Q)) = [];
    mean_Q(isinf(mean_Q)) = [];
end
if exist('gcFrac_1','var')
    % Fraction of nodes in giant component ignoring orphan nodes
%     mean_gcFrac_1 = sum(gcFrac_1)./sum(numMeasures100);
    gcFrac_1(isnan(gcFrac_1)) = 0;
    gcFrac_1(isinf(gcFrac_1)) = 0;
    for iD = 1:nM100
        mean_gcFrac_1(iD) = median(nonzeros(gcFrac_1(:,iD)));
    end
%     std_gcFrac_1 = sqrt(sum((gcFrac_1 - numMeasures100.*repmat(mean_gcFrac_1,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_gcFrac_1(isinf(std_gcFrac_1)) = 0;
%     std_gcFrac_1(isnan(std_gcFrac_1)) = 0;
%     std_gcFrac_1(isinf(mean_gcFrac_1)) = [];
    mean_gcFrac_1(isinf(mean_gcFrac_1)) = [];
    % Fraction of nodes in giant component ignoring orphan nodes
    mean_gcFrac_5 = sum(gcFrac_5)./sum(numMeasures100);
%     std_gcFrac_5 = sqrt(sum((gcFrac_5 - numMeasures100.*repmat(mean_gcFrac_5,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_gcFrac_5(isinf(std_gcFrac_5)) = 0;
%     std_gcFrac_5(isnan(std_gcFrac_5)) = 0;
%     std_gcFrac_5(isinf(mean_gcFrac_5)) = [];
    mean_gcFrac_5(isinf(mean_gcFrac_5)) = [];
    % Fraction of nodes in giant component ignoring orphan nodes
    mean_gcFrac_10 = sum(gcFrac_10)./sum(numMeasures100);
%     std_gcFrac_10 = sqrt(sum((gcFrac_10 - numMeasures100.*repmat(mean_gcFrac_10,optimax,1)).^2)./(sum(numMeasures100) - 1));
%     std_gcFrac_10(isinf(std_gcFrac_10)) = 0;
%     std_gcFrac_10(isnan(std_gcFrac_10)) = 0;
%     std_gcFrac_10(isinf(mean_gcFrac_10)) = [];
    mean_gcFrac_10(isinf(mean_gcFrac_10)) = [];
end
