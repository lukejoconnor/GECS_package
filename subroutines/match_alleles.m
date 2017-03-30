function [ phase ] = match_alleles( cell_array, allele1, allele2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

letters=[{'A'},{'T'},{'G'},{'C'}];
phase=cellfun(@(s,a1,a2)match(s,a1,a2),cell_array, allele1, allele2);


    function phase=match(s,a1,a2)
        if ~ismember(s,letters)
            phase=0;
        elseif ismember(s,[{a1},{a2}])
            phase=(-1)^(s==a1);
        else
            phase=(-1)^(comp(s)==a1);
        end
    end

    function let2=comp(let1)
        if let1=='A'
            let2='T';
        elseif let1=='T'
            let2='A';
        elseif let1=='C'
            let2='G';
        elseif let1=='G'
            let2='C';
        end
    end
end