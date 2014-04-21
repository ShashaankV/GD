%unpack the chromosome 22 correlation matrix if necessary
function unpack_chr22_SNPcorrelations(GDpath)
fo='chr22_SNPcorrelationmatrix'; %file name to save correlation matrix to 
                                    %if changed, confirm consistency with synthetic script 
%check if file exists, if it does then skip unpacking
tmp=dir(GDpath);
f={};
for i=1:length(tmp)
    b=tmp(i).name;
    f{i}=b(1:end-4); %strip *.mat extension
end
if ismember(fo,f)
    return 
else
a1=load('chr22_SNPcorrelations_1');
a2=load('chr22_SNPcorrelations_2');
Rlist=cat(1,a1.('R'),a2.('R')); %Rlist is the lower triangle (incl. diagonal), c1=row index, c2=column index, c3= correlation 
n=max(Rlist(:,1)); %get the number of parameters, SNPs
R=zeros(n); %initialize matrix
i=int32(Rlist(:,1)); %necessary to avoid MATLAB running into approximation issues for large indices
j=int32(Rlist(:,2));
n=int32(n);
ind=i+(j-1)*n;
R(ind)=Rlist(:,3);
R=R+tril(R,-1)';
save([GDpath '/' fo],'R','-v7.3')
end
end
