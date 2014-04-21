function [x,y]=y_additive(G,miss_value,h2,xs,xtype,xsign,varargin)
    if nargin>6
        xmaf=varargin{1};
    end
    [m,n]=size(G);
    A=zscore_sv(G,miss_value,'zero'); %this zscores (across subject rows) excluding missing genotypes and replaces missing with 0, 
                            %the third argument defines the value of the missing genotypes after normalization ('zero'=0,'random'=random value from the appropriate distribution)
                            %a fourth argument is optional, if 'T' then theoretical mean and standard deviation is used (only valid for Binomial G, ignored otherwise), default is empirical 
    
    %define shape of coefficients
    %if not hyperexponential then below defaults to Uniform nonzeros
    x=ones(n,1);
    if strcmp(xtype,'Hyperexponential1')
        x(:,1)=exp(-[1:n]/(.05*xs))+exp(-[1:n]/n); %make length constant xs
    elseif strcmp(xtype,'Hyperexponential2')
        x(:,1)=exp(-[1:n]/(.2*xs))+exp(-[1:n]/n); %make length constant xs
    end

    %assign signs
    if strcmp(xsign,'posneg')
        x(1:2:end,1)=-x(1:2:end,1); %flip every other sign
    elseif strcmp(xsign,'neg')
        x=-x; %flip all signs
    end
    
    %get rank order index in terms of maf if valid or randomly assign rank
    %order of indices
    if sum(sum(G(:,1)-floor(G(:,1))))==0 %test if elements in the first column of G are integers, if true assume Binomial matrix
        maf=calc_maf(G,miss_value); %calculate maf
        %assign location of nonzero coefficients
        if strcmp(xmaf,'random')
            ind=randperm(n);
        elseif strcmp(xmaf,'maf_low')
            %sort indices from high to low
            [i,ind]=sort(maf,'ascend');
        elseif strcmp(xmaf,'maf_high')
            [i,ind]=sort(maf,'descend');
        end
    else
        ind=randperm(n); %if not Binomial matrix then ignore maf criterion and randomly assign nonzero indices
    end
    
    %redefine x
    %1. set the first 1 to xs indices in the sorted index 'ind' of the new vector 'x' to the current 1
    %to xs values of 'x'...this reassigns the nonzero coefficients to these
    %locations in rank order according to the above
    x(ind(1:xs),1)=x(1:xs,1); 
    %2. filter the new 'x' through a mask to keep assignments as above but zero out what wasn't
    %moved...this also institutes the sparse criterion
    mask=zeros(n,1); %make all other zero
    mask(ind(1:xs),1)=1;
    x=x.*mask;
    
    %generate y, including scaling of x according to A and h2
    sigg=sqrt(h2);
    sige=sqrt(1-h2);
    sig=std(A*x);
    x=sigg*x/sig;
    z=sige*randn(m,1);
    y=A*x+z;
end
