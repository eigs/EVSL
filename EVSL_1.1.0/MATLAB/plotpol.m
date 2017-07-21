function [] = plotpol(mu, ab, gam, lmin, lmax)
dd = (lmax - lmin)/2;
cc = (lmax + lmin)/2;
ab0 = ab;
ab = (ab - cc)/ dd;
%% for visualizing the filter 
xi = [-1:0.01:1];
yi = jac2Pltd(mu, xi);
t = max(yi);
yi = yi/t;
%in = find((xi<ab(1))|(xi>ab(2)));
%tout = max(abs(yi(in)));
%mu = jac2Coefd(mbest, gam, 2);
%xi = [-1:0.01:1];
%yi = jac2Pltd(mu, xi);
%t = max(yi);
%yi = yi/t;
plot(xi,yi,[-1 1], [0 0],[gam gam],[0, 1],'linewidth',2)
hold on
plot([ab(1), ab(1)],[0, 1],'k--')
plot([ab(2), ab(2)],[0, 1],'k--')
end