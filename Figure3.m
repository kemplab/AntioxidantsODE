%% Sensitivity Analysis Code

% User inputs

% For ODE system inputs
% blap_ODE_system_proposal(endtimes, timestep, cell_density, protein_abund, blap_conc, p_change, x_change)
endtimes = 3600*2; % 2 hours in seconds
timestep = 1; % in seconds
cell_density = 1e9; % cells per Liter
blap_conc = 3;
p_change = ones(36,1);
x_change = ones(32,1);

%% Establish genes of interest and index
genes = ["GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", "PRDX1", "PRDX2", ...
    "PRDX3", "PRDX6", "CAT", "TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3", "GLRX", "GLRX2", ...
    "G6PD", "GLUD1", "GSR", "GSTP1", "POR", "NQO1", "SOD1", "SOD2", "SOD3", "AQP3", ...
    "AQP8", "AQP9", "GCLC"];
    
%% Get files
metadata_file = pwd() + "\\Code\\Data\\2021-06-08_scrnaseqdata.csv";
files = dir(pwd() + "\\Code\\Data\\HNSCC_Broad_scRNAseq");
%%
metadata = readtable(metadata_file);
%%
gene_index = ones(length(genes),1);
protein = readtable(files(3).folder + "\\" + metadata.cell(1) + ".csv");
for i = 1:length(genes)
    gene_index(i,1) = find(protein.Var1 == genes(i));
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
%%
warning_sims = csvread("20230112_failed_odes.csv");
%%
protein_vals(:,warning_sims) = [];
%%
protein_abund = mean(protein_vals, 2);

%% Run control sim and set control output variables

% establish number of parameters
num_params = 36;
% Sensitivities of all parameters

p_h2o2e_sensitivities = ones(num_params,2);
p_nadph_sensitivities = ones(num_params,2);
p_trx_sensitivities = ones(num_params,2);
p_gsh_sensitivities = ones(num_params,2);

y_ctrl = blap_ODE_system(endtimes, timestep, cell_density, protein_abund, blap_conc, p_change, x_change);
h2o2e_ctrl = y_ctrl(end,3);
nadph_ctrl = y_ctrl(end,24)/y_ctrl(end,25);
trx_ctrl = y_ctrl(end,15)/y_ctrl(end,16);
gsh_ctrl = y_ctrl(end,7)/y_ctrl(end,8);

%%
percentchanges = 10;
for percentchange = percentchanges

for param1 = 1:num_params
    %for kd_rate = 1:length(kd_rates)
    %disp(i)
    
    p_up = p_change;
    p_down = p_change;
    p_up(param1,1) = 1+percentchange/100;
    p_down(param1,1) = 1-percentchange/100;
    
    y_up =  blap_ODE_system(endtimes, timestep, cell_density, protein_abund, blap_conc, p_up, x_change);
    y_down =  blap_ODE_system(endtimes, timestep, cell_density, protein_abund, blap_conc, p_down, x_change);
    
    h2o2e_up = y_up(end,3);
    h2o2e_down = y_down(end,3);
    p_h2o2e_sensitivities(param1,1) = ((h2o2e_up-h2o2e_ctrl)/h2o2e_ctrl)/(percentchange/100);
    p_h2o2e_sensitivities(param1,2) = ((h2o2e_down-h2o2e_ctrl)/h2o2e_ctrl)/(percentchange/100);
    
    nadph_up = y_up(end,24)/y_up(end,25);
    nadph_down = y_down(end,24)/y_down(end,25);
    p_nadph_sensitivities(param1,1) = ((nadph_up-nadph_ctrl)/nadph_ctrl)/(percentchange/100);
    p_nadph_sensitivities(param1,2) = ((nadph_down-nadph_ctrl)/nadph_ctrl)/(percentchange/100);
    
    trx_up = y_up(end,15)/y_up(end,16);
    trx_down = y_down(end,15)/y_down(end,16);
    p_trx_sensitivities(param1,1) = ((trx_up-trx_ctrl)/trx_ctrl)/(percentchange/100);
    p_trx_sensitivities(param1,2) = ((trx_down-trx_ctrl)/trx_ctrl)/(percentchange/100);
    
    gsh_up = y_up(end,7)/y_up(end,8);
    gsh_down = y_down(end,7)/y_down(end,8);
    p_gsh_sensitivities(param1,1) = ((gsh_up-gsh_ctrl)/gsh_ctrl)/(percentchange/100);
    p_gsh_sensitivities(param1,2) = ((gsh_down-gsh_ctrl)/gsh_ctrl)/(percentchange/100);
end
end
%%
save_file = 1;
t_sens = 3600*2;
p_filename = '2023-01-12_p_fig3.png';

% Use below for single parameter sensitivity

p_names = {'k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12', ...
    'k13','k14','k15','k16','k17','k18','k19','k20','k21','k22','k23', ...
    'k24','k25','k26','k27','k28','k29','k30','k31','k32','k33','k34', 'k35', 'k36'};

str_pctchange="10";
for sens = 1:4
    % H2O2 Sensitivity
    figure('Position', [10 10 800 1500])
    if sens == 1
        Objective_low_sum = p_h2o2e_sensitivities(:,2);
        Objective_high_sum = p_h2o2e_sensitivities(:,1);
        filename = "p_h2o2e_sens_final_"+date+"_"+str_pctchange+".png";
    
    % NADPH sens
    elseif sens == 2
        Objective_low_sum = p_nadph_sensitivities(:,2);
        Objective_high_sum = p_nadph_sensitivities(:,1);
        filename = "p_nadph_sens_final_"+date+"_"+str_pctchange+".png";
        
    % GSH sens
    elseif sens == 3
        Objective_low_sum = p_gsh_sensitivities(:,2);
        Objective_high_sum = p_gsh_sensitivities(:,1);
        filename = "p_gsh_sens_final_"+date+"_"+str_pctchange+".png";
    
    % Trx sens
    elseif sens == 4
        Objective_low_sum = p_trx_sensitivities(:,2);
        Objective_high_sum = p_trx_sensitivities(:,1);
        filename = "p_Trx_sens_final_"+date+"_"+str_pctchange+".png";
    end
    % Sort the values based on the lower change
    % Sort the higher values and the names arrays
    %    using the same indices
    if sens == 1
        [~,ind]=sort(Objective_high_sum,'ascend');
        names_Objective = p_names(ind);
    end
    Objective_low_sum = Objective_low_sum(ind);
    Objective_high_sum = Objective_high_sum(ind);
    % Create a figure and plot the low and high horizontally
    h = barh(Objective_high_sum);
    hold on
    xmax=max([max(Objective_low_sum),max(Objective_high_sum),-min(Objective_low_sum),-min(Objective_high_sum)]);
    xlim([-1.025*xmax 1.025*xmax])
    barh(Objective_low_sum,'r')
    if sens == 1
        title("Intracellular $[H_{2}O_{2}]$ Sensitivities",'interpreter','latex',"FontSize", 24);
        xlabel("$\frac{\Delta{}[H_{2}O_{2}]}{\Delta{}k}$",'interpreter','latex','FontSize',24);
    elseif sens == 2
        title("Intracellular $\frac{[NADPH]}{[NADP^+]}$ Sensitivities",'interpreter','latex',"FontSize", 24);
        xlabel("$\frac{\Delta{}\frac{[NADPH]}{[NADP^+]}}{\Delta{}k}$",'interpreter','latex','FontSize',24);
    elseif sens == 3
        title("Intracellular $\frac{[GSH]}{[GSSG]}$ Sensitivities",'interpreter','latex',"FontSize", 24);
        xlabel("$\frac{\Delta{}\frac{[GSH]}{[GSSG]}}{\Delta{}k}$",'interpreter','latex','FontSize',24);
    elseif sens == 4
        title("Intracellular $\frac{[rTrx]}{[oTrx]}$ Sensitivities",'interpreter','latex',"FontSize", 24);
        xlabel("$\frac{\Delta{}\frac{[rTrx]}{[oTrx]}}{\Delta{}k}$",'interpreter','latex','FontSize',24);
    end
    
    set(gca,'Ytick',[1:length(names_Objective)],'YTickLabel',[1:length(names_Objective)])
    set(gca,'yticklabel',names_Objective)

    % get the current tick labels
    ticklabels = get(gca,'YTickLabel');
    % prepend a color for each tick label
    ticklabels_new = cell(size(ticklabels));

    
    for i = [3, 4, 5, 13, 20, 24, 25, 26, 35]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{red} ' ticklabels{ordered_i}];
    end
    for i = [8, 9, 10, 11, 12, 21, 27, 28]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{blue} ' ticklabels{ordered_i}];
    end
    for i = [7, 22]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color[rgb]{0, 0.8, 0} ' ticklabels{ordered_i}];
    end
    for i = [14, 15, 16, 17, 18, 19]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color[rgb]{1,0.6,0.2} ' ticklabels{ordered_i}];
    end
    for i = [29, 30, 31, 32, 33, 34, 36]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{black} ' ticklabels{ordered_i}];
    end
    for i = [1, 2, 6, 23]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{magenta} ' ticklabels{ordered_i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_new, "FontSize", 20);
    
    %xlabel('colors for legend')
    l1 = plot([NaN,NaN], 'r');
    l2 = plot([NaN,NaN], 'b');
    l3 = plot([NaN,NaN], 'color', '#00CC00');
    l4 = plot([NaN,NaN], 'color', '#FF9933');
    l5 = plot([NaN,NaN], 'k');
    l6 = plot([NaN,NaN], 'm');
    [~,hObj]=legend([l1, l2, l3, l4, l5, l6], {'GSH', 'Prx', ...
        'NADPH', 'Protein Thiol', 'Drug Metabolism', 'Other'},...
        'Position', [.18,.35,.3,.1], "FontSize", 20);           % return the handles array
    legend boxoff
    hL=findobj(hObj,'type','line');  % get the lines, not text
    set(hL,'linewidth',1.5)
    
    if save_file
        saveas(gcf,filename)
    end
end
