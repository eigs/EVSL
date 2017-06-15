function [yi] = jac2Pltd (mu, xi)
%% function [yi] = jac2Pltd (mu, gam, damping, xi)
%% computes the values of pn(xi) given the
%% expansion coefficients muA
%%-------------------- compute p(xi)
n = size(xi,1) ;
yi = zeros(n,1);
vkm1 = zeros(n,1);
vk   = ones(n,1);
for k=0:length(mu)-1;
    %%-------------------- accumulation of vector.
    yi = yi + mu(k+1)*vk;
    scal = 2;
    if (k == 0) , scal = 1;, end;
    vkp1 = scal .*xi .* vk - vkm1;
    vkm1 = vk;
    vk = vkp1;
end
