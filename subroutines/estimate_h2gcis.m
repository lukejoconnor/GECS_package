%% RUNS GECOSCORE REGRESSION AND ESTIMATE DENOMINATOR
% 8/31/16.

%% Load + process data
no_cis_snps=6900;

% Subsample 1
p=1;
load([qtl_path,num2str(p),'.mat'])
allele1=a1;allele2=a2;
SNPs=uniqueSNPs;
beta_eqtl(beta_eqtl.^2>10)=0;
beta_prod=beta_eqtl;

% Subsample 2
p=2;
load([qtl_path,num2str(p),'.mat'])
[intersection,idx1,idx2]=intersect(SNPs,uniqueSNPs);
no_snps=length(intersection);
if size(beta_eqtl,1)~=size(beta_prod,1)
    error('gene mismatch between subsamples')
end
beta_eqtl(beta_eqtl.^2>10)=0;
beta_prod=beta_prod(:,idx1).*beta_eqtl(:,idx2);

% LDSC data
load(ldsc_path)
[~,idxldsc,idx0]=intersect(rsid,intersection,'stable');
no_snps=length(idx0);
LDscores=LDscores(idxldsc,1);
beta_prod=beta_prod(:,idx0);

%% Estimate cis h2g

ldsc_mat=(beta_prod~=0)*spdiags(LDscores,0,no_snps,no_snps);
h2gcis_estimate=no_cis_snps*sum(ldsc_mat.*beta_prod,2)./sum(ldsc_mat.^2,2);
h2gcis=mean(h2gcis_estimate(h2gcis_estimate==h2gcis_estimate));

% Jackknife blocks: exclude every SNP pair s.t. either SNP is in excluded block
if no_blocks>0   
    exclude=sparse(no_snps,no_blocks,no_snps);
    blocklength=floor(no_snps/no_blocks);
    for kk=1:no_blocks
        exclude((1+blocklength*(kk-1)):(blocklength*kk),kk) = true;
    end
    exclude((1+blocklength*(kk-1)):(no_snps),no_blocks) = true;
end

for jk=1:no_blocks
    temp=no_cis_snps*sum(ldsc_mat(:,~exclude(:,jk)).*beta_prod(:,~exclude(:,jk)),2)./sum(ldsc_mat(:,~exclude(:,jk)).^2,2);
    h2gcis_jk_estimate(jk)=mean(temp(temp==temp));
end

