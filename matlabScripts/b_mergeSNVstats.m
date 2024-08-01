%% merges null and obs SNVstats per chromosome

%% Executed in parallel using available threads
disp('Starting ... mergeSNVstats');
parfor cChr = 1:22
    fprintf('# chromosome: chr%d\n', nChr);
    
    % Set the observed and null input files per chromosome
    obs_fname = ['../SNVstats/' cohortName '.obs.chr' num2str(cChr) '.mat'];
    null_fname = ['../SNVstats/' cohortName '.null.chr' num2str(cChr) '.mat'];
    
    % Read in the observed data
    obs_data = load(obs_fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell','N_samp','N_snv');
    snv_ids = obs_data.snv_ids;
    snv_refs = obs_data.snv_refs;
    snv_alts = obs_data.snv_alts;
    samp_ids = obs_data.samp_ids;
    sampXsnv_cell = obs_data.sampXsnv_cell;
    N_samp = null_data.N_samp;
    N_snv = null_data.N_snv;

    % Read in the null data
    null_data = load(null_fname,'snv_ids','snv_refs','snv_alts','samp_ids','sampXsnv_cell','N_samp','N_snv');
    snv_ids2 = null_data.snv_ids;
    snv_refs2 = null_data.snv_refs;
    snv_alts2 = null_data.snv_alts;
    samp_ids2 = null_data.samp_ids;
    sampXsnv_cell2 = null_data.sampXsnv_cell;
    N_samp2 = null_data.N_samp;
    N_snv2 = null_data.N_snv;
    
    % Merge the observed and null
    %snv_ids = [snv_ids snv_ids2];
    snv_ids = {snv_ids{1:N_snv} snv_ids2{1:N_snv2}};
    snv_refs = [snv_refs snv_refs2];
    snv_alts = [snv_alts snv_alts2];
    %samp_ids = [samp_ids, strcat(samp_ids2, '-null')];
    samp_ids = {samp_ids{1:N_samp} samp_ids2{1:N_samp2}};
    %sampXsnv_cell = [sampXsnv_cell, cellfun(@(x) x + N_snv, sampXsnv_cell2, 'UniformOutput', false)];
    sampXsnv_cell = {sampXsnv_cell{1:N_samp} sampXsnv_cell2{1:N_samp2}};

    for i = 1:N_samp2
        samp_ids{N_samp+i} = [samp_ids{N_samp+i} '-null'];
        sampXsnv_cell{N_samp+i} = sampXsnv_cell{N_samp+i} + N_snv;
    end

    cCmp = {snv_ids{N_snv+1:length(snv_ids)}};
    for i = 1:N_snv
        vec = strcmp(snv_ids{i},cCmp);
        if sum(vec) > 0
            display([num2str(i) '/' num2str(N_snv)]);
            idx = find(vec==1);
            idx1 = i; idx2 = idx + N_snv;
            snv_ids = {snv_ids{1:idx2-1} snv_ids{idx2+1:length(snv_ids)}};
            snv_refs = [snv_refs(1:idx2-1) snv_refs(idx2+1:end)];
            snv_alts = [snv_alts(1:idx2-1) snv_alts(idx2+1:end)];
            for j = 1:N_samp2
                vec2 = sampXsnv_cell{N_samp+j};
                vec2(vec2==idx2) = idx1;
                vec2(vec2>idx2) = vec2(vec2>idx2) - 1;
                sampXsnv_cell{N_samp+j} = vec2;
            end
            cCmp = {snv_ids{N_snv+1:length(snv_ids)}};
        end
    end
    
    % update stats
    N_snv = numel(snv_ids);
    N_samp = numel(samp_ids);
    fprintf('# snvs: %d\n', N_snv);
    fprintf('# samps: %d\n', N_samp);
    snv_shared = zeros(1, N_snv);
    for i = 1:N_samp
        snv_shared(sampXsnv_cell{i}) = snv_shared(sampXsnv_cell{i}) + 1;
    end
    h_shared = histcounts(snv_shared, 1:max(snv_shared)+1);
    
    fname2 = strrep(obs_fname, 'obs.chr', 'obs.null.merged.chr');
    s = struct(...
        'snv_ids', {snv_ids}, ...
        'snv_refs', snv_refs, ...
        'snv_alts', snv_alts, ...
        'samp_ids', {samp_ids}, ...
        'sampXsnv_cell', {sampXsnv_cell}, ...
        'N_samp', N_samp, ...
        'N_snv', N_snv, ...
        'snv_shared', snv_shared, ...
        'h_shared', h_shared ...
    );
    save(fname2, "-fromstruct", s);
end
