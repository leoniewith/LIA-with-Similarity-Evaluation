function [H,p_value,Z] = Mann_Kendall(V)
% p_value <= 0.05 有显著性
% H = 1 有趋势
% H =0 没有趋势
% Z > 0 上升趋势
% Z < 0 下降趋势
alpha = 0.05;
V = reshape(V,length(V),1);
alpha = alpha/2; %
n = length(V);
i = 0; 
j = 0; 
S = 0;
for i = 1:n-1
    for j = i+1:n
        S = S + sign(V(j) - V(i));
    end
end
VarS = (n*(n-1)*(2*n+5))/18;
StdS = sqrt(VarS);
%%%% Note: ties are not considered
if S >= 0
    Z = ((S-1)/StdS)*(S~=0);
else
    Z = (S+1)/StdS;
end
p_value = 2*(1-normcdf(abs(Z),0,1)); %% Two-tailed test
pz = norminv(1-alpha,0,1);
H = abs(Z)>pz; %%
return
end