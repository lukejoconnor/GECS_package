function process_qtls( path, frq_path, save_path )
% process_qtls processes eQTL summary statistics files in FastQTL format
%  path (string): path to file containing eQTL files ([path, '*.qtl']).
%  frq_path (string): path to file containing allele count data in
%  PLINK .hwe format, generated using --hardy. If desired to use only a
%  subset of SNPs in downstream analyses, exclude some SNPs from this file
%  (e.g. using option --extract *.in)
%  save_path (string): save path ('*.mat')

if nargin~=3
    error('Wrong number of inputs')
end


% Files in directory that match
temp=dir([path,'*.qtl']);
if isempty(temp)
    error('Bad eQTL file path')
end
for i=1:length(temp) tissues{i}=temp(i).name; end

% Load eQTL data
counter=0; %Keep track of number of entries loaded so far
genecounter=0; %Number of gene-tissue pairs loaded so far
for celltypeindex=1:length(tissues)
    
    fprintf('%d ',celltypeindex)
    
    % Import eQTL data
    [betat,genest,SNPst ] = import_eQTL_sumstats( [path,tissues{celltypeindex}] );%2nd argument: no. entries to load
    no_entries=length(SNPst);
    [genest,~,geneidxt]=unique(genest);
    
    % Allocating memory
    if celltypeindex==1
        SNPs=zeros(round(no_entries*length(tissues)*1.1),1);
        rowidx=SNPs;
        beta=SNPs;
        genes=cell(round(max(geneidxt)*1.1),1);
        rowtissues=zeros(round(max(geneidxt)*1.1),1);
    end
    
    % Recording data for current tissue
    SNPs(counter+1:counter+no_entries)=SNPst;
    rowidx(counter+1:counter+no_entries)=genecounter+geneidxt;
    beta(counter+1:counter+no_entries)=betat;
    genes(genecounter+1:genecounter+max(geneidxt))=genest;
    rowtissues(genecounter+1:genecounter+max(geneidxt))=celltypeindex;
    
    counter=counter+no_entries;
    genecounter=genecounter+max(geneidxt);
    
end
clear *t

fprintf('\n')

SNPs=SNPs(1:counter);
rowidx=rowidx(1:counter);
beta=beta(1:counter);
genes=genes(1:genecounter);
rowtissues=rowtissues(1:genecounter);

[uniqueSNPs, ~, colidx]=unique(SNPs);
no_snps=length(uniqueSNPs);
no_genes=genecounter;


% Matrix of eqtl effects.
beta_eqtl=sparse(rowidx,colidx,beta,no_genes,no_snps);
clear rowidx colidx beta

% Allele phase and allele freq data
file=fopen(frq_path);fgets(file);
data=textscan(file, '%*d rs%d %*s %s %s %d/%d/%d %*f %f %*f\n');

[uniqueSNPs,idx_qtl,idx_ref]=intersect(uniqueSNPs,data{1},'stable');

beta_eqtl=beta_eqtl(:,idx_qtl); 
no_snps=length(idx_qtl);

% Char for ref allele and alt allele    
a1=data{2}(idx_ref);
a2=data{3}(idx_ref);

% Allele frequency
heterozygosity=double(data{7}(idx_ref));

% Convert beta_eqtl from per-allele to per-s.d.
beta_eqtl=beta_eqtl*diag(sparse(sqrt(heterozygosity)));


% Saving
save(save_path,'beta_eqtl','a1','a2','no_snps','no_genes','uniqueSNPs',...
    'genes','rowtissues','heterozygosity','-v7.3')


end

