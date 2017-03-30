% Loads LD data computed from each subsample of the eQTL cohort

% QTL data
load([qtl_path,num2str(0),'.mat'])
beta_eqtl_orig=beta_eqtl;

[no_genes,no_snps_orig]=size(beta_eqtl_orig);


% Discard outliers in beta_eqtl (only seem to exist for subsample?)
beta_eqtl_orig(beta_eqtl_orig.^2>1)=0;

intersection{1}=uniqueSNPs;

% Load LD data for each subsample
for p=1:2
    file=fopen([ld_path,num2str(p),'.ld']);fgets(file);
    LDdata=textscan(file,'%d %d rs%d %*d %*d rs%d %f\n');
    if ~feof(file)
        warning('Failed to load LD data from whole file')
    end
    fclose(file);
    
    % Build sparse matrix
    [LDdata{3},mapback1,idxld1]=unique(LDdata{3});
    [LDdata{4},~,idxld2]=unique(LDdata{4});
    ldcs4_orig{p}=sparse(idxld1,idxld2,LDdata{5});
    
    % Make rows and columns match
    [RSIDld{p},idxld1,idxld2]=intersect(LDdata{3},LDdata{4},'stable');
    ldcs4_orig{p}=ldcs4_orig{p}(idxld1,idxld2);
    
    % Reorder rows and columns (for jackknifing)
    chr=LDdata{1}(mapback1(idxld1));
    pos=LDdata{2}(mapback1(idxld1));
    [~,reorder]=sortrows([chr pos]);
    ldcs4_orig{p}=ldcs4_orig{p}(reorder,reorder);
    RSIDld{p}=RSIDld{p}(reorder);
    
    % Match SNPs with LD data
    [intersection{p+1}, idxld{p}, idx3{p}]=intersect(RSIDld{p},intersection{p},'stable');
    ldcs4_orig{p}=ldcs4_orig{p}(idxld{p},idxld{p});
    [iii,jjj,entries]=find(ldcs4_orig{p});
    RR_orig{p}=sparse(iii,jjj,entries);
    clear iii jjj entries 
    
    ldcs4_orig{p}=ldcs4_orig{p}.^2;
    
    ldcs4_orig{p}=max(ldcs4_orig{p},ldcs4_orig{p}') - 1/n_subsample*(max(ldcs4_orig{p},ldcs4_orig{p}')~=0);
    if p==2
        ldsc_inv=1./max(1,sum(ldcs4_orig{p})');
    end
    ldcs4_orig{p}=triu((ldcs4_orig{p})^2 );
    
    clear LDdata
end

beta_eqtl_orig=beta_eqtl_orig(:,idx3{1}(idx3{2}));
ldcs4_orig{1}=ldcs4_orig{1}(idx3{2},idx3{2});

RR_max_orig=max(RR_orig{1}(idx3{2},idx3{2}).^2,RR_orig{2}.^2);
RR_max_orig=max(RR_max_orig,RR_max_orig');

RR_orig=triu(RR_orig{1}(idx3{2},idx3{2})+RR_orig{2})/2;% Unsigned LD matrix
heterozygosity_orig=heterozygosity;
ldsc_inv_orig=ldsc_inv;
a1_orig=a1(idx3{1}(idx3{2}));

% Build shared LD matrix
ldcs4_orig=max(triu(ldcs4_orig{1}+ldcs4_orig{1}'),triu(ldcs4_orig{2}+ldcs4_orig{2}'));

clear temp*
intersection_orig=intersection{3};

% For jackknifing
idxkeep=idxld{1}(idx3{2});
idxkeep_orig=idxkeep;


