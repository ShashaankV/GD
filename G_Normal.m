function G=G_Normal(m,n,miss,miss_value,R)
G=randn(m,n);
U=cholcov(R);
if size(U,1)*size(U,2)==n^2
    G=G*U;
%     C2=corrcoef(G);
%     c2=C2(triu(ones(n),1)~=0);
%     R=R(triu(ones(n),1)~=0);
%     plot(R,c2,'.')
else
    fprintf('U is not square. Correlation matrix of G will be identity')
end


[m,n]=size(G);
for i=1:n
    for j=1:m
        if rand>1-miss %missing genotype probability (1=no missing)
            G(j,i)=miss_value;
        end
    end
end

end