function G=G_Binomial(m,n,miss,miss_value,maf)
G=zeros(m,n);
for i=1:n
    x=binornd(2,maf(i),[m,1]);
    G(:,i)=x;
end
maf=reshape(maf,[n,1]);

[m,n]=size(G);
for i=1:n
    for j=1:m
        if rand>1-miss %missing genotype probability (1=no missing)
            G(j,i)=miss_value;
        end
    end
end

end