%% Settings and setup
format compact

initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

%% Read transcriptomics data
% genes were mapped from gene name to Entrez in Python with the myGene module
% df = tdfread('Data/microarray_data_with_entrez_genes.csv'); % Microarray
df = tdfread('./Data/RNAseq_data_with_entrez_genes_Cell_AVERAGE(missing=1).txt'); % RNA-seq

% store the cell line columns individually
Average = df.Average;

datasetnames = {'Average'};
datasets = {Average};

%% Calculate RAS scores for all reactions
R3 = readCbModel('Models/medium_only/Recon3DModel_301_patched_M1.mat'); % load one of our patched medium models

RASmatrix = zeros(length(R3.rxns),length(datasets));
for i = 1:length(datasets)
    data = datasets{i};
    transcripts = data;
    
    % loop over all reactions to define the RAS
    for j=1:length(R3.rxns)
        if(strcmp(R3.rules(j),'')) % no genes linked to this reaction
            Coeff = NaN;
        else
            rules = R3.rules{j};
            Coeff = getScore(rules, transcripts); 
        end
        RASmatrix(j,i) = Coeff; 
    end
end

% How many reactions will get a restricted bound?
length(find(~isnan(RASmatrix(:,1))))

% export the data for analysis in Python
save('Data/RAS_scores.mat','RASmatrix','-mat')

%% Plot maximal biomass flux vs. NNR factor on all 18 models
% which NNR factors to try?
scales = logspace(3,-1,50); % logarithmically from 1000 to 0.1
final_NNR_index = 13; % was 13
final_NNR = scales(final_NNR_index);

% store results in a matrix
% columns for the NNR factor, rows for the 12 models
% values will be the maximal biomass flux
results = zeros(18,length(scales));

% loop over the two datasets
model_count = 1;
for j = 1:length(datasets)
    % get RAS scores for this cell line
    cell_line = datasetnames{j};
    RAScolumn = RASmatrix(:,j); % PLC, Huh7

    % loop over the 6 medium models
    for medium = {'1','2','3','4','5','6'}
        % load model
        boundedModel_orig = readCbModel(['Models/medium_only/Recon3DModel_301_patched_M' char(medium) '.mat']);

        % loop over NNR factor
        NNR_count = 1;
        for NNR = scales
            boundedModel = boundedModel_orig;
            maxScore = 1000;

            % loop over all rxns and adjust bounds
            for i=1:length(boundedModel.rxns)
                rxn_val = RAScolumn(i);

                if isnan(rxn_val) % no genes linked to reaction
                    % check if bounds are 1000 (this excludes specifically set
                    % exchange reactions
                    if boundedModel.lb(i) == -1000
                        boundedModel.lb(i) = -1*maxScore;
                    end
                    if boundedModel.ub(i) == 1000
                        boundedModel.ub(i) = maxScore;
                    end
                else % nonzero value
                    if boundedModel.lb(i) == -1000
                        boundedModel.lb(i) = -1*NNR*rxn_val;
                    end
                    if boundedModel.ub(i) == 1000
                        boundedModel.ub(i) = NNR*rxn_val;
                    end
                end
            end

            % store maximal biomass flux
            sol = optimizeCbModel(boundedModel,'max');
            results(model_count, NNR_count) = sol.f;

            % save models for the final NNR factor
            if NNR == final_NNR 
                writeCbModel(boundedModel, 'mat', ['./Models/Recon3DModel_301_patched_M' char(medium), '_', cell_line, '_NNR=', num2str(NNR),'.mat']);
            end
            
            NNR_count = NNR_count + 1;
            
        end
        
        model_count = model_count + 1;
    end
end

% plot results
figure('pos',[10 10 1000 400],'DefaultAxesFontSize',16)
semilogx(scales,results(1:6,:),'--','LineWidth', 2) % HL60
% hold on;
% ax = gca;
% ax.ColorOrderIndex = 1;
% semilogx(scales,results(7:12,:),':','LineWidth', 2) % Huh7
% hold on;
% ax.ColorOrderIndex = 1;
% semilogx(scales,results(13:18,:),'-.','LineWidth', 2) % Patient
xlim([-0.1,1000])
xlabel('NNS')
ylabel('Maximal biomass flux')
legend('HL60 M1','HL60 M2','HL60 M3','HL60 M4','HL60 M5','HL60 M6')
line([final_NNR, final_NNR],[0,0.00001],'Color','black','LineStyle','--','HandleVisibility','off')
hold off;

fig = gcf; fig.PaperPositionMode = 'auto'; % on-screen size
export_fig('Figures/biomass_vs_NNR.png','-png','-r200','-p0.01')

% biomass bar plot
figure('pos',[10 10 1000 400],'DefaultAxesFontSize',16)
bar(results(1:6,final_NNR_index),'k')
xticklabels({'HL60 M1','HL60 M2','HL60 M3','HL60 M4','HL60 M5','HL60 M6'})
xtickangle(45)
ylabel('Maximal biomass flux')
export_fig('Figures/biomass_at_final_NNR.png','-png','-r200','-p0.01')


%% Write all models to SBML
% for cell_line = {'PLC','Huh7'}
%     for medium = {'1','2','3','4','5','6'}
%         modelName = ['../Models/Recon3DModel_301_patched_M' char(medium), '_', char(cell_line), '_NNR=', num2str(scales(final_NNR_index)),'.mat'];
%         m = readCbModel(modelName);
%     writeCbModel(m, 'sbml', modelName);
%     end
% end