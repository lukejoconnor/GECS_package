function [ ] = import_LDscores( path_to_ldsc, save_path, path_to_snpcounts )
%Processes LDSC file to generate LD scores.

chr=[];
rsid=[];
pos=[];
maf=[];
LDscores=[];

for i=1:22
    file=fopen([path_to_ldsc,num2str(i),'.l2.ldscore']);fgets(file);
    data=textscan(file, '%d rs%d %d %*f %f %f\n');
    if ~feof(file)
        error('failed to parse data fully')
    end
    chr=[chr; data{1}];
    rsid=[rsid; data{2}];
    pos=[pos; data{3}];
    maf=[maf; data{4}];
    LDscores=[LDscores; data{5}];
    
    fclose(file);
end

if exist('path_to_snpcounts')
    file=fopen(path_to_snpcounts);
    data=textscan(file,'%d\n');
    total_no_snps=sum(data{1});
    save(save_path,'chr','rsid','pos','maf','LDscores','total_no_snps');
else
    save(save_path,'chr','rsid','pos','maf','LDscores');
end



end

