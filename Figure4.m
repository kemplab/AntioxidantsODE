timestep = 1; % in seconds
endtimes = 3600*2; % 2 hours in seconds

% Basal
p_change = ones(36,1);
%p_change(1,1) = 0;
x_change = ones(32,1);
blap_exists = 1; 
cell_density = 1e9; % cells per Liter

%% Establish genes of interest
genes = ["GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", ...
    "PRDX1", "PRDX2", "PRDX3", "PRDX6", "CAT", "TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3", ...
    "GLRX", "GLRX2", "G6PD", "GLUD1", "GSR", "GSTP1", "POR", "NQO1", "SOD1", "SOD2", "SOD3", ...
    "AQP3", "AQP8", "AQP9", "GCLC"];
    
%% Get files
metadata_file = pwd() + "\\Code\\Data\\2021-06-08_scrnaseqdata.csv";
files = dir(pwd() + "\\Code\\Data\\HNSCC_Broad_scRNAseq");
%%
metadata = readtable(metadata_file);
%%
gene_index = [];
protein = readtable(files(3).folder + "\\" + metadata.cell(1) + ".csv");
for i = 1:length(genes)
    gene_index = [gene_index, find(protein.Var1 == genes(i))];
end
%%
protein_vals = zeros(33,height(metadata));
%%
for i = 1:height(metadata)
    protein = readtable(files(3).folder + "\\" + metadata.cell(i) + ".csv");
    protein_vals(:,i) = table2array(protein(gene_index,2));
    
end
% Convert from PPM to micromolar
protein_vals = protein_vals.*5.29e-4;

%% Run control sim
protein_averages = mean(protein_vals, 2);

%%
protein_vals(34,:) = metadata.Tumor;
protein_vals(35,:) = metadata.non_cancerIndex;


%%
results_blap = zeros(32,length(protein_vals));
%%
%figure
tic
warning_sims = [];
for sim = 1:length(protein_vals)
    if mod(sim,50) == 0
        disp(sim)
        toc
    end
    %results_ctrl(:,sim) = y_ctrl(end,2:end);
    % reset lastwarn
    lastwarn('', '');
    % call the ode solver that fails sometimes
    y_blap = blap_ODE_system(endtimes, timestep, cell_density, protein_vals(:,sim), 3, p_change, x_change);
    % lastwarn now has info
    [warnMsg, warnId] = lastwarn();
    % only keep simulations that do not throw error and save others as
    % warning sims
    if(isempty(warnId))
        results_blap(:,sim) = y_blap(end,2:end);
    else
        warning_sims = [warning_sims, sim];
    end
    
end

%%

full_results = [results_blap;
    protein_vals];

csvwrite("20230112_ode_results.csv",full_results)
csvwrite("20230112_failed_odes.csv",warning_sims)
            
