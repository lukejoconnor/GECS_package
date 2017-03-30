% Wrapper script for GEcoscore regression after some intermediary files
% have been computed. Requires:
%   - Signed GWAS summary statistics in .sumstats format. Order of columns
%   must be SNPs, Z scores, allele 1, allele 2, N, *. Path to file
%   containing GWAS sumstats is phenotype_path
%   - Processed eQTL summary statistics generated using process_qtls.m,
%   both for whole sample and for 2 subsamples. Filenames should be
%   [qtl_path,'.0.mat'], [qtl_path,'.1.mat'], [qtl_path,'.2.mat'] for
%   whole-sample and subsample summary stats respectively
%   - LD scores from reference panel (run import_LDscores on text file)
%   - In-sample LD data for each eQTL subsample. This should be generated
%   using e.g.: plink2 --bfile path --r --ld-window 999999 --ld-window-kb 
%   $windowsize_in_kb --remove $samples_to_exclude --out
%   ld_path_$subsample_number

% Path to GWAS summary stats ([phenotype_path,'*.sumstats'])
phenotype_path='.../';

% Path to LDscore data
ldsc_path='... .mat';

% Path to processed eQTL data
qtl_path='...';

% Path to LD data
ld_path='...';
load_ld_matfile=false; % Option to save LD data as a .mat file for efficiency

% Where to save - will be appended with phenotype number and '.mat'
path_to_save='...';

% Results can be generated for different parameter settings
quantiles=[0 .25 .5];% Proportion of SNP pairs to retain
weight_regression=true(1,3); % Weight regression to improve power?
heterozygosity_lb=.05*ones(1,3);% Lower bound on heterozygosity for SNPs
heterozygosity_ub=ones(1,3);% Upper bound on heterozygosity for SNPs
correct_regression=true(1,3);% Condition on LD (in addition to pruning)?

no_blocks=50; % Number of jackknife blocks

total_no_snps=9.25*10^6;% Number of 1kg SNPs, for h2gtrait estimation
avg_no_cis_snps=6900;% Average number of cis SNPs, for h2gcis estimation

% Names of GWAS datasets
temp=dir([phenotype_path,'*.sumstats']);
for i=1:length(temp) phenotypes{i}=temp(i).name; end
no_pheno=length(phenotypes);

% Load eQTL and LD data
if load_ld_matfile
    load(ld_path);
else
    load_LDdata;
end
fprintf('Finished loading LD data\n')

% Estimate h2gcis using cross-val
estimate_h2gcis
fprintf('Estimated h2gcis: %f\n',h2gcis)

% Loop over phenotypes
pheno_array=1:no_pheno;
for jj=pheno_array
    pheno_path=[phenotype_path, phenotypes{jj}];
    disp(phenotypes{jj});
    save_path=[path_to_save,num2str(jj),'.mat'];

    % Run GEcoscore regression.
    GEcoscore_regression;
    
end