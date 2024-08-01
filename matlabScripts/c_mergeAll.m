%% merges all chromosome level SNVstats into one file

%% Executed in parallel using available threads
disp('Starting ... mergeAll');
fname_all = ['../SNVstats/' cohortName '.obs.null.merged.mat'];
snv_ids_all = cell(0,0);
snv_refs_all = [];
snv_alts_all = [];
samp_ids_all = cell(0,0);
sampXsnv_cell_all = cell(0,0);
N_samp_all = 0;
N_snv_all = 0;

parfor cChr = 1:22
    fprintf('# chromosome: chr%d\n', cChr);
    
    fname = strrep(fname_all,'.mat',['.chr' num2str(cChr) '.mat']);
    
    % Read in the per chromosome data
    chrom_data = load(fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell','N_samp','N_snv');
    snv_ids = chrom_data.snv_ids;
    snv_refs = chrom_data.snv_refs;
    snv_alts = chrom_data.snv_alts;
    samp_ids = chrom_data.samp_ids;
    sampXsnv_cell = chrom_data.sampXsnv_cell;
    N_samp = chrom_data.N_samp;
    N_snv = chrom_data.N_snv;

    % Create temporary variables to store results
    temp_snv_ids_all = snv_ids;
    temp_snv_refs_all = snv_refs;
    temp_snv_alts_all = snv_alts;
    temp_samp_ids_all = samp_ids;
    temp_sampXsnv_cell_all = sampXsnv_cell;

    N_snv0 = length(snv_ids_all);
    N_samp0 = length(samp_ids_all);
    for i = 1:N_samp
        temp_sampXsnv_cell_all{i} = temp_sampXsnv_cell_all{i} + N_snv0;
    end
    
    % Combine data into temporary variables
    snv_ids_all = [snv_ids_all, temp_snv_ids_all];
    snv_refs_all = [snv_refs_all, temp_snv_refs_all];
    snv_alts_all = [snv_alts_all, temp_snv_alts_all];

    if isempty(samp_ids_all)
        samp_ids_all = temp_samp_ids_all;
        sampXsnv_cell_all = temp_sampXsnv_cell_all;
    else
        for i = 1:data.N_samp
            idx = find(strcmp(temp_samp_ids_all{i}, samp_ids_all));
            if isempty(idx)
                samp_ids_all = [samp_ids_all, temp_samp_ids_all{i}];
                sampXsnv_cell_all = [sampXsnv_cell_all, temp_sampXsnv_cell_all{i}];
            else
                sampXsnv_cell_all{idx} = [sampXsnv_cell_all{idx}, temp_sampXsnv_cell_all{i}];
            end
        end
    end
end

snv_ids = snv_ids_all;
snv_refs = snv_refs_all;
snv_alts = snv_alts_all;
samp_ids = samp_ids_all;
sampXsnv_cell = sampXsnv_cell_all;

% update stats
N_snv = length(snv_ids);
N_samp = length(samp_ids);
fprintf('# snvs: %d\n', N_snv);
fprintf('# samps: %d\n', N_samp);
snv_shared = zeros(1, N_snv);
for i = 1:N_samp
    snv_shared(sampXsnv_cell{i}) = snv_shared(sampXsnv_cell{i}) + 1;
end
h_shared = histcounts(snv_shared, 1:max(snv_shared)+1);

save(fname_all,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell',...
    'N_samp','N_snv','snv_shared','h_shared');
