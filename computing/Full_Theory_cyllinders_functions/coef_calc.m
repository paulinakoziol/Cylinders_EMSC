function [A, B, C, D, V, W] = coef_calc(l_max, xi, eta, n_hat, zeta, ...
    bessel_eta, bessel_xi, hankel_xi, bessel_diff_eta, bessel_diff_xi, hankel_diff_xi)
% Calculations of intermediate parameters based on: 
% Applied Spectroscopy 2019, Vol. 73(8) 859–869

l = (0:l_max);
l = reshape(l, 1, 1, size(l,2));
l = l.*ones(size(xi));

A = 1i*xi.*(xi.*bessel_diff_eta.*bessel_xi-eta.*bessel_eta.*bessel_diff_xi);

B = xi.*(n_hat.*n_hat.*xi.*bessel_diff_eta.*bessel_xi-eta.*bessel_eta.*bessel_diff_xi);

C = l.*cos(zeta).*eta.*bessel_eta.*bessel_xi.*((xi.*xi)./(eta.*eta)-1);

D = l.*cos(zeta).*eta.*bessel_eta.*hankel_xi.*((xi.*xi)./(eta.*eta)-1);

V = xi.*(n_hat.*n_hat.*xi.*bessel_diff_eta.*hankel_xi-eta.*bessel_eta.*hankel_diff_xi);

W = 1i.*xi.*(eta.*bessel_eta.*hankel_diff_xi-xi.*bessel_diff_eta.*hankel_xi);

end

