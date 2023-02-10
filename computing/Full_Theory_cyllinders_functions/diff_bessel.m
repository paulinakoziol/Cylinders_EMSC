function [bessel_eta, bessel_xi, hankel_xi, bessel_diff_eta, bessel_diff_xi, ...
    hankel_diff_xi] = diff_bessel(l_max, xi, eta)
% This function calculates bessel and hankel functions along with their derivatives, 
% for l from 0 to l_max range and for arguments xi and eta 

l = (0:l_max);
l = reshape(l, 1, 1, size(l,2));
l = l.*ones(size(xi));

% Bessel and Hankel functions
bessel_eta = besselj(l, eta.*ones(1,1,l_max+1));
bessel_xi = besselj(l, xi.*ones(1,1,l_max+1));
hankel_xi = bessel_xi+1i*bessely(l, xi.*ones(1, 1, l_max+1));

%Bessel derivatives
bessel_diff_eta = differentiate_num(bessel_eta, eta);
bessel_diff_xi = differentiate_num(bessel_xi, xi);
%Hankel derivative
hankel_diff_xi = differentiate_num(hankel_xi, xi);

end

