function [ SNPs, Z, a1, a2, N  ] = import_sumstats( path_to_file )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

file=fopen(path_to_file);
[~,temp]=fscanf(file,'%[SNP] %[A1] %[A2] %[N] %[CHISQ] %[Z]');
if temp~=6
    error('Column headings not in expected order, or columns missing')
end

data=textscan(file,'rs%d %c %c %f %*f %f\n');
SNPs=data{1};
a1=cellstr(data{2});
a2=cellstr(data{3});
N=data{4};
Z=data{5};

if ~feof(file)
    warning('Did not import all SNP entries from summary statistics file')
end

end

