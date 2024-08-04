%% merges all chromosome level SNVstats into one file

%% Executed in parallel using available threads
disp('Starting ... mergeAll');

% Output filename
fname_all = ['../SNVstats/' cohortName '.obs.null.merged.mat'];

% Set variables to hold final merged data
snv_ids_all = cell(1, 22);
snv_refs_all = cell(1, 22);
snv_alts_all = cell(1, 22);
samp_ids_all = cell(1, 22);
sampXsnv_cell_all = cell(1, 22);
N_samp_all = zeros(1, 22);
N_snv_all = zeros(1, 22);

parfor cChr = 1:22
    fprintf('# chromosome: chr%d\n', cChr);
    
    % Per chromosome input filename
    fname = strrep(fname_all,'.mat',['.chr' num2str(cChr) '.mat']);
    
    % Read in the per chromosome data
    chrom_data = load(fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell','N_samp','N_snv');
    snv_ids_all{cChr} = chrom_data.snv_ids;
    snv_refs_all{cChr} = chrom_data.snv_refs;
    snv_alts_all{cChr} = chrom_data.snv_alts;
    samp_ids_all{cChr} = chrom_data.samp_ids;
    sampXsnv_cell_all{cChr} = chrom_data.sampXsnv_cell;
    N_samp_all(cChr) = chrom_data.N_samp;
    N_snv_all(cChr) = chrom_data.N_snv;

end

% Collapse the data per chromosome into single variable
snv_ids = [snv_ids_all{:}];
snv_refs = [snv_refs_all{:}];
snv_alts = [snv_alts_all{:}];
samp_ids = unique([samp_ids_all{:}]);
sampXsnv_cell = cell(1, length(samp_ids));

% Need to adjust per chromosome sampXsnv_cell indices 
N_snv_cumsum = [0, cumsum(N_snv_all)];
for cChr = 1:22
    for i = 1:length(samp_ids_all{cChr})
        samp_idx = find(strcmp(samp_ids_all{cChr}{i}, samp_ids));
        sampXsnv_cell{samp_idx} = [sampXsnv_cell{samp_idx}, sampXsnv_cell_all{cChr}{i} + N_snv_cumsum(cChr)];
    end
end


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

disp('Completed ... mergeAll');
