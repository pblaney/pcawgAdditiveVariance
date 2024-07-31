%%% additive variance of passengers in MM

cohortName = 'Multiple-Myeloma';

cd matlabScripts

%%% call matlab pipeline
a_makeKeys
b_mergeSNVstats
c_mergeAll
d_makeOrderedKey
e_makeMACHfiles
f_call_gcta
f_summarize_results

cd ..
