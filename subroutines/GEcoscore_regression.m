%% Runs GE coscore regression
% 11/26/16

%% Load + process data

% GWAS data
[ SNPph, Zph, allele1, allele2, Nph  ] = import_sumstats( pheno_path );


% Match eQTL and GWAS data with LD data
[intersection, idx1, idx2]=intersect(intersection_orig,SNPph,'stable');

no_snps=length(idx1);

Nph=double(Nph(idx2));
Zph=Zph(idx2)./sqrt(Nph);%Phenotype sumstats

% Alt/ref allele phasing needs to be consistent (ok if everything is flipped)
allele1=allele1(idx2);allele2=allele2(idx2);
phase=match_alleles(a1_orig(idx1),allele1,allele2);

beta_eqtl=beta_eqtl_orig(:,idx1)*diag(sparse(phase));
RR=diag(sparse(phase))*RR_orig(idx1,idx1)*diag(sparse(phase));
RR_max=RR_max_orig(idx1,idx1);
ldcs4=ldcs4_orig(idx1,idx1);

heterozygosity=heterozygosity_orig(idx1);
clear phase_inv phase


fprintf('matched %d SNP IDs\n',no_snps);

fprintf('%d out of %d genes have eQTLs\n',full(sum(any(beta_eqtl'~=0))),no_genes)

clear SNPph temp RSIDs SNPph RSIDld phase coscores

p=1;
gecs{p}=triu(beta_eqtl'*beta_eqtl);
gecs{p}=gecs{p}-diag(diag(gecs{p}));% Remove diagonal entries
idxkeep=idxkeep_orig;
a1=a1_orig;

% Subsample and its complement (phase of subsample doesn't matter)

for p=1:2
    load([qtl_path,num2str(p),'.mat'],'beta_eqtl','uniqueSNPs','a1')
    
    % Discard outliers in beta_eqtl (why do these outliers exist?!)
    beta_eqtl(beta_eqtl.^2>1)=0;
    
    % Gene mismatch with subsample
    if (size(beta_eqtl,1)~=no_genes)
        warning('Gene count mismatch between subsample + whole sample')
    end
    
    % Handle SNP mismatch with subsample
    [intersection,tidx2,tidx1]=intersect(uniqueSNPs,intersection,'stable');
    
    a1=a1(tidx2);% for phasing
    ldcs4=ldcs4(tidx1,tidx1);
    RR=RR(tidx1,tidx1);
    RR_max=RR_max(tidx1,tidx1);
    
    ldsc_inv=ldsc_inv_orig(tidx1);
    idxkeep=idxkeep(tidx1);
    allele1=allele1(tidx1);allele2=allele2(tidx1);
    Zph=Zph(tidx1);
    for pp=1:p
        gecs{pp}=gecs{pp}(tidx1,tidx1);
    end
    no_snps=length(intersection);
    
    
    % Subsample phase
    phase=match_alleles(a1,allele1,allele2);
    beta_eqtl=beta_eqtl(:,tidx2)*spdiags(phase,0,(no_snps),(no_snps));
    
    gecs{p+1}=triu(beta_eqtl'*beta_eqtl);
    gecs{p+1}=gecs{p+1}-diag(diag(gecs{p+1}));% Remove diagonal entries
    
end

clear beta_eqtl

% Get SNP pairs with nonzero gecs in all 3 samples
nzcs=gecs{1}~=0;
for p=2:3
    nzcs=nzcs.*(gecs{p}~=0);
end
nzcs=nzcs~=0;
[iii,jjj]=find(nzcs);

regression_weights_nonuniform=ldsc_inv(iii).*ldsc_inv(jjj);


for p=1:3
    coscores{p}=gecs{p}(nzcs);
    gecs{p}=[];
end

ldcs4=triu(ldcs4+ldcs4');%%%
RR=triu(RR+RR');%%%
RR_max=triu(RR_max+RR_max');%%%

[~,temp]=sort(full(ldcs4(nzcs)));
ranking1(temp)=1:nnz(nzcs);
[~,temp]=sort(full(RR_max(nzcs)));
ranking2(temp)=1:nnz(nzcs);
clear distances

if exist('LD_case')
    if LD_case==1 % Shared-LD only
        distances=-ranking1;
    elseif LD_case==2 % LD only
        distances=-ranking2;
    else % Combined
        [~,temp]=sort(max(ranking1,ranking2));
        distances(temp)=-(1:nnz(nzcs));
    end
else
    [~,temp]=sort(max(ranking1,ranking2));
    distances(temp)=-(1:nnz(nzcs));
end
distances=distances';
clear ranking1 ranking2

RR=full(RR(nzcs));


bandsizes=quantile(distances,quantiles);

% Minimum and maximum heterozygosity for each SNP pair
min_het=min(heterozygosity(iii),heterozygosity(jjj));
max_het=max(heterozygosity(iii),heterozygosity(jjj));


% Jackknife blocks: exclude every SNP pair s.t. either SNP is in excluded block
exclude=cell(no_blocks,1);
blocklength=floor(max(idxkeep)/no_blocks);
if no_blocks>0
    for kk=1:no_blocks
        exclude{kk} = (find((idxkeep(iii)>=1+blocklength*(kk-1)).*(idxkeep(iii)<=blocklength*kk)...
            +(idxkeep(jjj)>=1+blocklength*(kk-1)).*(idxkeep(jjj)<=blocklength*kk)));
    end
    exclude{no_blocks} = (find((idxkeep(iii)>=1+blocklength*(kk-1)).*(idxkeep(iii)<=max(idxkeep))...
        +(idxkeep(jjj)>=1+blocklength*(kk-1)).*(idxkeep(jjj)<=max(idxkeep))));
end

%% Estimate regression denominator

% Loop over LD thresholds
for bb=1:numel(bandsizes)
    keep=(min_het>=heterozygosity_lb(bb)).*(max_het<=heterozygosity_ub(bb));
    keep= 1==(keep.*(distances>=bandsizes(bb)));
    
    if weight_regression(bb)
        regression_weights=ldsc_inv(iii).*ldsc_inv(jjj);
    else
        regression_weights=ones(length(iii),1);
    end
    
    % Estimate for complement of each jackknife block
    for jk=1:no_blocks
        keepjk=false(length(keep),1);
        keepjk(exclude{jk})=keep(exclude{jk});
        if correct_regression(bb)
            for p=1:3
                correction_constant(p)=weighted_regression(RR(keepjk),coscores{p}(keepjk),regression_weights(keepjk));
            end
            block_estimate(jk,bb)=mean(...
                (coscores{2}(keepjk)-RR(keepjk)*correction_constant(2)).*...
                (coscores{3}(keepjk)-RR(keepjk)*correction_constant(3)).*...
                regression_weights(keepjk)) / mean(regression_weights(keepjk));
        else
            block_estimate(jk,bb)=mean(coscores{2}(keepjk).*coscores{3}(keepjk).*regression_weights(keepjk))/mean(regression_weights(keepjk));
        end
    end
    block_estimate(block_estimate~=block_estimate)=0;% region may have 0 SNPs if excluded from summary stats (eg MHC)
    block_weights=cellfun(@(x)sum(regression_weights(x)),exclude);
    
    % Point estimate
    if correct_regression(bb)
        for p=1:3
            correction_constant(p)=weighted_regression(RR(keep),coscores{p}(keep),regression_weights(keep));
        end
        denominator_estimate(bb)=mean(...
            (coscores{2}(keep)-RR(keep)*correction_constant(2)).*...
            (coscores{3}(keep)-RR(keep)*correction_constant(3)).*...
            regression_weights(keep)) / mean(regression_weights(keep));
    else
        denominator_estimate(bb)=mean(coscores{2}(keep).*coscores{3}(keep).*...
            regression_weights(keep)) / mean(regression_weights(keep));%sum(weights.*block_estimate(:,bb))/sum(weights);
    end
    
   
    % Jackknife estimates: take each block back out
    for jk=1:no_blocks
        jk_denominator_estimate(jk,bb)=denominator_estimate(bb)-block_estimate(jk,bb)*block_weights(jk)/sum(block_weights);
    end
end

coscores=coscores{1};
clear block_estimate beta_eqtl_array gecs

jk_denominator_estimate=mean(jk_denominator_estimate,3);
denominator_estimate=mean(denominator_estimate,1);
denominator_pct_err=std(jk_denominator_estimate)*sqrt(no_blocks+1)./denominator_estimate;
if any(denominator_pct_err>1/3)
    warning('High percent error in denominator, leading to possibly unstable estimates. Possibly your effective eQTL sample size is too small.:')
    disp((denominator_pct_err))
end
if any(denominator_pct_err<0
    warning('Negative denominator estimates. Possibly your effective eQTL sample size is too small.')
end
%% Run GEcoscore regression

% Products of phenotypic effect sizes
Zcoscores=Zph(iii).*Zph(jjj);

% GEcoscore regression
for bb=1:numel(bandsizes)
    
    if weight_regression(bb)
        regression_weights=ldsc_inv(iii).*ldsc_inv(jjj);
    else
        regression_weights=ones(length(iii),1);
    end
    
    compute_GEcoscore_numerator
    
    jk_error(bb)=sig2alpha_std;
    unnormalized_estimates(bb)=sig2alpha_est;
    unnormalized_jk_estimates(bb,:)=jk_est;
    
end

% Unnormalized estimates: not normalized by total variance ratio for phenotype : gexp
estimates=unnormalized_estimates./denominator_estimate;
jk_estimates=unnormalized_jk_estimates'./jk_denominator_estimate;

%% Estimate trait heritability for denominator

% Load and match SNPs
load(ldsc_path,'rsid','LDscores')
[~,idxldsc,idx0]=intersect(rsid,intersection,'stable');
LDscores=LDscores(idxldsc,1);
ldsc_inv=ldsc_inv(idx0);
idxkeep=idxkeep(idx0);% For jackknifing
Zph=Zph(idx0);

temp=weighted_regression([LDscores ones(length(LDscores),1)],Zph.^2,ldsc_inv);
h2g_pheno=temp(1)*total_no_snps;%%%

% Jackknifing
keepidxinv=zeros(blocklength*no_blocks,1);
keepidxinv(idxkeep)=1:length(idxkeep);
for jk=1:no_blocks
    keep=nonzeros(keepidxinv([1:blocklength*(jk-1), blocklength*jk:blocklength*no_blocks]));
    temp=weighted_regression([LDscores(keep) ones(length(keep),1)],Zph(keep).^2,ldsc_inv(keep));
    h2g_jk_pheno(jk)=temp(1)*total_no_snps;
end

estimates=estimates*h2gcis/h2g_pheno;
jk_estimates=diag(h2gcis_jk_estimate./h2g_jk_pheno)*jk_estimates;
jk_error=std(jk_estimates)*sqrt(no_blocks+1);

%% Saving
save([save_path,num2str(jj),'.mat'],...
    'unnormalized_jk_estimates','unnormalized_estimates',...
    'jk_denominator_estimate','jk_estimates','estimates','jk_error',...
    'h2g_jk_pheno','h2g_pheno','h2gcis','h2gcis_jk_estimate')

beep
%%
%save(['concat_estimates_radius=',num2str(cis_radius/1000),'kb.mat'],'estimates','jk_error','phenotypes');
