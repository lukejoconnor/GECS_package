% Runs GEcoscore "regression" without the regression part
% Regression denominator (i.e. var(gecs)) must be estimated separately

%% Get point estimates

include_pair=(min_het>=heterozygosity_lb(bb)).*(max_het<=heterozygosity_ub(bb));
include_pair= 1==(include_pair.*(distances>=bandsizes(bb)));


% Point estimate
if correct_regression(bb)
    correction_constant=weighted_regression(RR(include_pair),Zcoscores(include_pair),regression_weights(include_pair));
    
    temp=mean(coscores(include_pair).*...
        (Zcoscores(include_pair)-RR(include_pair)*correction_constant) .*...
        regression_weights(include_pair)) / mean(regression_weights(include_pair));
    [ sig2alpha_est ]=temp*no_genes;
else
    temp = mean(coscores(include_pair).*Zcoscores(include_pair).*regression_weights(include_pair))/mean(regression_weights(include_pair));
    [ sig2alpha_est ]=temp*no_genes;
end





%% Jackknifing
if no_blocks>0
    jk_est=zeros(no_blocks,1);
    include_pair0=include_pair;
    for ii=1:no_blocks
        include_pair=include_pair0;
        include_pair(exclude{ii})=false;
        
        if correct_regression(bb)
            correction_constant=weighted_regression(RR(include_pair),Zcoscores(include_pair),regression_weights(include_pair));
            
            temp=mean(coscores(include_pair).*...
                (Zcoscores(include_pair)-RR(include_pair)*correction_constant) .*...
                regression_weights(include_pair)) / mean(regression_weights(keep));
        else
            temp = mean(coscores(include_pair).*Zcoscores(include_pair).*regression_weights(include_pair))/mean(regression_weights(include_pair));
        end
        
        [ jk_est(ii) ]=temp*no_genes;
        
    end
    jk_estimate=sum(jk_est/no_blocks);
    sig2alpha_std=sqrt((no_blocks-1)/no_blocks*sum((jk_est-jk_estimate).^2));
    
else
    sig2alpha_std=0;
end

%sig2alpha_std_covar=sqrt((no_blocks-1)/no_blocks*sum((jk_est2-jk_estimate).^2));

