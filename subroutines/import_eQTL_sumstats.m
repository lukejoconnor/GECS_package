function [ beta_per_allele, genes, snps, pval] = import_eQTL_sumstats( path )
%import_gtex_downsamp loads eQTLs in the fastQTL file format
%   Called by process_qtls

file=fopen(path);

data=textscan(file,'%s rs%d %*d %f %f\n');

if ~feof(file)
    error(['Failed to parse whole file. Filename: ',path])
end

genes=data{1};
snps=data{2};
pval=data{3};
beta_per_allele=data{4};

end

