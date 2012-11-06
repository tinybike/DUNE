% ppi.m is an implementation of the DUNE (DUplication & NEofunctional-
% ization) model of eukaryotic protein-protein interaction network 
% evolution.  Please see http://ppi.tinybike.net for project details.  To 
% run, this script requires several supporting files:
%
% ppi_import.m (data import, requires adjacency matrix as input file)
% pointmutation.m (neofunctionalization/assimilation events)
% timedepstats.m (updates dynamical properties at specified intervals)
% trackstats.m (updates degree and betweenness "tracking" values)
% 
% In addition to these, the MatlabBGL package:
% http://www.mathworks.com/matlabcentral/fileexchange/10922
% and Louvain modularity algorithm:
% https://sites.google.com/a/brain-connectivity-toolbox.net/bct/
% are also required.
%
% The following files can be used to calculate various statistics and plot
% several static and dynamic quantities obtained from the simulation, but
% are not required to run ppi.m:
%
% std_dyn.m (calculates dynamical statistics)
% std_stat.m (calculates static end-state statistics)
% std_tracking.m (calculates degree and betweenness "tracking" statistics)
% ppi_dyn.m (prepares dynamical quantities for plotting)
% RT_save_figures_bvk.m (betweenness vs degree)
% RT_save_figures_nvk.m (nearest-neighbor degree vs degree)
% RT_save_figures_Cvk.m (clustering coefficients vs degree)
% RT_save_figures_E2.m (second-largest eigenvalue evolution)
% RT_save_figures_px.m (closeness distribution)
% RT_save_figures_plambda.m (walk matrix eigenvalue distribution)
% RT_save_figures_pb.m (betweenness distribution)
% RT_save_figures_pk.m (degree distribution)
% RT_save_figures_gcc.m (global clustering coefficient evolution)
% RT_save_figures_Q.m (modularity coefficient evolution)
% RT_save_figures_d.m (diameter evolution)
% RT_save_figures_f1.m (largest component evolution)
%
% Many PPI data files are stored as edge lists.  To convert HitPredict
% edge lists to adjacency matrices, we used hitpredict_A.m.  This script
% is picky about the edge list format.  Columns must be as follows:
% Protein 1 ID | Protein 2 ID | Concatenated ID's | Unique ID's only
% (We used Microsoft Excel to format our edge lists this way.)
% 
% For more flexible edge list manipulation (including tools for converting
% edge lists to adjacency matrices), I recommend visiting Gergana Bounova's
% website:
% http://www.mit.edu/~gerganaa/downloads.html
%
% Please direct any questions/comments to jack@tinybike.net.
%
% (Note, in the paper, the parameters shown below use different names:
% a -> c, mu -> p_mutation, d -> p_duplicate, and phi -> p_silence)

best_RMSE = Inf;
savecount = 0;

clearvars -except data_*

% Program variables
optimax = 50;           % Number of replicate simulations
numBinsResults = 100;   % Number of bins for degree sequence
numBins = 100;          % Number of bins for other sequences
frags = 20;             % Number of fragmentation data points
max_neighbors = 20;     % Maximum number of neighbors for assimilation
tmax = 100000;          % Maximum number of time steps
N_0 = 2;                % Initial number of proteins

% NOTE: In the paper, the parameters use different names
% (a -> c, mu -> p_mutation, d -> p_duplicate, and phi -> p_silence)

% Fruit fly
c = 0.54608;
p_mutation = 0.00058781;
p_duplicate = 0.0014;
p_silence = 0.86632;
 
% Yeast
c = 0.68979;
p_mutation = 0.0007861;
p_duplicate = 0.01;
p_silence = 0.5546;

% Human
c = 0.727274;
p_mutation = 0.000762;
p_duplicate = 0.003700;
p_silence = 0.651732;

p_A = p_duplicate/(p_mutation + p_duplicate);

% Import experimental data
ppi_import

fprintf('%i (%i) nodes, %i (%i) edges in L1\n',data_N,nnz(sum(data_A)),data_edges,full(sum(sum(data_A)) + trace(data_A))/2);

% Define bins to ensure bin centers are the same for all sequences
avgBinsResults = linspace(1,100,numBinsResults);
avgBins = linspace(0,1,numBins);
avgBinsB = linspace(0,max(data_b)*1.5,numBins);
avgBinsCls = linspace(0,max(data_cls)*2,numBins);
avgBinsEigs = linspace(-1,1,numBinsResults);
avgBinsNV = linspace(0,max(data_NV)*1.5,numBins);

% Define empty arrays for histograms
total_y = zeros(1,numBinsResults);
total_cy = zeros(1,numBins);
total_by = zeros(1,numBins);
total_clsy = zeros(1,numBins);
total_ey = zeros(1,numBinsResults);
total_nvy = zeros(1,numBins);
total_knny = zeros(1,numBinsResults);

% Define empty arrays for scatter plots
maxUniqueK = 1000;
total_mean_C = zeros(maxUniqueK,1);
total_mean_b = zeros(maxUniqueK,1);
total_mean_knn = zeros(maxUniqueK,1);
total_mean_cls = zeros(maxUniqueK,1);
tau = zeros(optimax,tmax);

% Define empty cell arrays to store raw results
k_cell = cell(optimax,1);
C_cell = cell(optimax,1);
b_cell = cell(optimax,1);
x_cell = cell(optimax,1);
A_cell = cell(optimax,1);

% Outer loop: replicate simulations
for opt_m = 1:optimax
    
    % Initialize adjacency matrix A with two connected nodes
    N = N_0;
    A = sparse(zeros(N));
    A(1,2) = 1;
    A(2,1) = 1;
    B = sparse(zeros(N));
    timeDep = 1;
    mut_accum = zeros(N,1);
	
	% Inner loop: time steps
    for t = 1:tmax
        
        % Real-time simulation (Gillespie)
        % tau: time step duration
        g_rand = rand(1);
        tau(opt_m,t) = 1/N/(p_mutation + p_duplicate)*log(1/rand(1));
        
        [linker,linker] = max(rand(1,N));
        
        % Point mutation (neofunctionalization)
        if g_rand >= p_A
            A = pointmutation(A,linker,c,max_neighbors);
        % Duplication/divergence
        else
            A(N+1,:) = A(linker,:);
            A(:,N+1) = cat(1,A(N+1,:)',0);
            duplicateEdges = find(A(linker,:));
            numCopy = length(duplicateEdges);
            if ~isempty(duplicateEdges)
                for di = 1:numCopy
                    if rand(1) < p_silence
                        pickRand = rand(1);
                        if pickRand > 0.5
                            A(linker,duplicateEdges(di)) = 0;
                            A(duplicateEdges(di),linker) = 0;
                        else
                            A(N+1,duplicateEdges(di)) = 0;
                            A(duplicateEdges(di),N+1) = 0;
                        end
                    end
                end
            end
            N = N + 1;
        end
		
        simEdges = (sum(sum(A)) + trace(A))/2;
        totalNumEdges(opt_m,t) = simEdges;
        totalNumNodes(opt_m,t) = N;
        avgNumEdges(opt_m,t) = totalNumEdges(opt_m,t)/N;
		
		% Calculate gene losses
		orphan(opt_m,t) = length(find(sum(A) == 0));
		orphanfrac(opt_m,t) = orphan(opt_m,t)/(length(A) - orphan(opt_m,t));
        
		% Tracking: uncomment to do individual-protein betweenness and
        % degree tracking
        trackstats;
        
        % Evolutionary trajectory calculations (diameter, modularity, etc.)
        if sum(sum(A)) > 0
            timedepstats;
        end
        
        % If simulated network edges >= experimental network edges, then
        % stop simulation and calculate statistics
        if simEdges >= (sum(sum(data_A)) + trace(data_A))/2
            break;
        end
                
    end
    
	% Calculate end-state statistics for comparison to data
    final_snapshot;
    
end