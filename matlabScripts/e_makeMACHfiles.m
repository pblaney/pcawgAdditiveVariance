%% writes files necessary for GCTA to gctaFiles folder (separate inputs for each funseq threshold)
%% inputs are from keys and bedfiles folders
%% additionally, an intermediate file is saved in the machMats folder, containing the genotype matrices saved as a matlab structure

disp('Starting ... makeMACHfiles');
input_tag = cohortName;
outFold = '../gctaFiles/';

bedFold = '../bedFiles/';
keyFold = '../keys/';
machFoldSNV = '../machMats/';
mergeFold = '../SNVstats/';

load([keyFold input_tag '.orderedKey.mat']);
load([mergeFold input_tag '.obs.null.merged.mat']);

nSNV = length(ordKey_fsq);
nSamp = N_samp;
filt = (snv_shared>1);
filt_idx = find(filt);
nFilt = sum(filt);
if ~(nFilt==length(ordKey_missing))
    display('filter test failed');
end

machMatsSNV = zeros(nSamp,nSNV,6);
machMat0SNV = zeros(nSamp,nSNV);
machMatsSNV_noncd = zeros(nSamp,nSNV,6);
machMat0SNV_noncd = zeros(nSamp,nSNV);
machMatsSNV_nonprm = zeros(nSamp,nSNV,6);
machMat0SNV_nonprm = zeros(nSamp,nSNV);
machMatsSNV_drv = zeros(nSamp,nSNV,6);
machMat0SNV_drv = zeros(nSamp,nSNV);
snv_chr_all = zeros(1,nSNV);
rej = 0;

for cSNV = 1:nSNV
    
    if mod(cSNV,1000)==0
        display([num2str(cSNV) '/' num2str(nSNV)]);
    end
    
    cIdx = filt_idx(cSNV);
    
    sampIdxs = [];
    for cSamp = 1:nSamp
        if ismember(cIdx,sampXsnv_cell{cSamp})
            sampIdxs = [sampIdxs cSamp];
        end
    end
    
    snv_chr = snv_ids{cIdx};
    clnIdx = find(snv_chr==':',1);
    snv_chr = snv_chr(4:clnIdx-1);
    snv_chr = str2num(snv_chr);
    if isempty(snv_chr)
        snv_chr = 23;
    end
    if snv_chr < 1
        snv_chr = 23;
    end
    snv_chr_all(cSNV) = snv_chr;
    
    snv_cd = ordKey_cd(cSNV);
    snv_drv = ordKey_drv(cSNV);
    snv_prm = ordKey_prm(cSNV);
    snv_fsq = ordKey_fsq(cSNV);
    lim = floor(snv_fsq);
    
    if lim>0
        if snv_drv==1
            machMatsSNV_drv(sampIdxs,cSNV,1:lim) = 1;
        else
            machMatsSNV(sampIdxs,cSNV,1:lim) = 1;
            if snv_cd==1
                machMatsSNV_noncd(sampIdxs,cSNV,1:lim) = 1;
                if snv_prm==0
                    machMatsSNV_nonprm(sampIdxs,cSNV,1:lim) = 1;
                end
            end
        end
    else
        if snv_drv==1
            machMat0SNV_drv(sampIdxs,cSNV) = 1;
        else
            machMat0SNV(sampIdxs,cSNV) = 1;
            if snv_cd==1
                machMat0SNV_noncd(sampIdxs,cSNV) = 1;
                if snv_prm==0
                    machMat0SNV_nonprm(sampIdxs,cSNV) = 1;
                end
            end
        end
    end
end

for i = 1:length(samp_ids)
    if ~isempty(strfind(samp_ids{i},'null'))
        samp_ids{i} = strrep(samp_ids{i},'-null','_null');
    end
end

singSamps = zeros(length(samp_ids),1);
phenVec = zeros(length(samp_ids),1);
for i = 1:length(samp_ids)
    if ~isempty(strfind(samp_ids{i},'null'))
        vec = strcmp(samp_ids{i}(1:end-5),samp_ids);
        if sum(vec)==0
            singSamps(i) = 1;
        end
        continue;
    else
        phenVec(i) = 1;
        vec = strcmp([samp_ids{i} '_null'],samp_ids);
        if sum(vec)==0
            singSamps(i) = 1;
        end
    end
end
phenVec = phenVec(singSamps==0);

fname = [machFoldSNV input_tag '.mat'];
save(fname,'machMatsSNV','machMat0SNV','machMatsSNV_noncd','machMat0SNV_noncd',...
    'machMatsSNV_nonprm','machMat0SNV_nonprm','machMatsSNV_drv','machMat0SNV_drv',...
    'samp_ids','rej','singSamps','phenVec','snv_chr_all','-v7.3');

% write out gcta files
output_tag = [outFold input_tag];
cMachMats = machMatsSNV;
cMachMat0 = machMat0SNV + machMatsSNV(:,:,1);

% Parallel execution of creating .dose / .info / .phen GCTA files
parfor cFsq = 0:6
    if cFsq==0
        cMat = cMachMat0;
    else
        cMat = cMachMats(:,:,cFsq);
    end

    fprintf('# FunSeq score threshold: %d\n', cFsq);
    cMat = cMat(singSamps==0,:);
    cBinMat = (cMat>0);
    colTots = sum(cBinMat,1);
    colInds = find((colTots>1)&(snv_chr_all<=22));
    cMat = cMat(:,colInds);
    nGene = size(cMat,2);
    N_samp = size(cMat,1);
    
    % Enhance the performance by writing the completed loop to the file instead of per line
    % Dose file
    fprintf('# Generating fsq%d.dose file\n', cFsq);
    dose_fid = fopen([output_tag '.fsq' num2str(cFsq) '.dose'],'w');
    dose_lines = strings(N_samp, 1);
    for i = 1:N_samp
        % Use of the array function to better execute the loop within a loop
        dose_lines(i) = sprintf('samp%d ALLELE_DOSE %s\n', i, strjoin(arrayfun(@(x) num2str(x, '%.4f'), cMat(i, :), 'UniformOutput', false), ' '));
    end
    fprintf(dose_fid, '%s', dose_lines);
    fclose(dose_fid);
    
    % Info file
    fprintf('# Generating fsq%d.info file\n', cFsq);
    info_fid = fopen([output_tag '.fsq' num2str(cFsq) '.info'],'w');
    fprintf(info_fid,'SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n');
    info_lines = strings(nGene, 1);
    info_maf = 0.5;
    info_freq = 0.5;
    for i = 1:nGene
        info_lines(i) = sprintf('snp%d\tA\tT\t%.4f\t%.4f\t1.0000\t1.0000\n', i, info_freq, info_maf);
    end
    fprintf(info_fid,'%s',info_lines);
    fclose(info_fid);
    
    % Phenotype file
    fprintf('# Generating fsq%d.phen file\n', cFsq);
    phen_fid = fopen([output_tag '.fsq' num2str(cFsq) '.phen'],'w');
    nonSings = find(singSamps==0);
    phen_lines = strings(N_samp, 1);
    for i = 1:N_samp
        ii = nonSings(i);
        phen_lines(i) = sprintf('samp%d\tsamp%d\t%d\n', i, i, ~contains(samp_ids{ii},'null'));
    end
    fprintf(phen_fid,'%s',phen_lines);
    fclose(phen_fid);
end

disp('Completed ... makeMACHfiles');
