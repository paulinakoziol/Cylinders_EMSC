function [Q] = full_theory_cylinders(WN, r, n_p, zeta, orientation)
% This function calculates scattering efficiency for samples with
% cylindrical domains using full Mie type Scattering Theory for cylinders
%
% Arguments description: 
% - WN - wavenumber vector
% - r - cyliner radius 
% - n_p - particles refractive index (real part) vector
% - zeta - ratiation incidence angle 
% - orientation - cylinder orientation (degrees) with respect to polarization 
%
% Calculations of Scattering Efficiency for samples with cylindrical
% domains is based on full Mie type scattering theory provided in:
% Applied Spectroscopy 2019, Vol. 73(8) 859–869

r = gpuArray(r);
n_p = gpuArray(n_p);

n_h = 1; 
WN = WN.*100; %conversion from cm^-1 to m^-1

x = 2*pi*n_h.*WN.*r;
x_hat = 2*pi.*WN.*r;
n_hat = n_p./n_h; %n_hat is defined as n_p/n_h but n_h is considered to be equal 1

xi = x_hat.*sin(zeta);
eta = x_hat.*sqrt(n_hat.*n_hat-cos(zeta).*cos(zeta));
l_max = x+4.05*(x).^(1/3)+2;
l_max = ceil(l_max); 
l_max = max(l_max, [], 2);
l_max = max(l_max);

% Calculations of bessel functions and their derivatives
[bessel_eta, bessel_xi, hankel_xi, bessel_diff_eta, bessel_diff_xi, ...
    hankel_diff_xi] = diff_bessel(l_max, xi,eta);

% Calculations of intermediate parameters
[A, B, C, D, V, W] = coef_calc(l_max, xi, eta, n_hat, zeta, bessel_eta, ...
    bessel_xi, hankel_xi, bessel_diff_eta, bessel_diff_xi, hankel_diff_xi);

% Scattering efficiency
Q = scattering_efficiency(x, A, B, C, D, V, W, orientation);
Q = gather(Q);

% figure;
% plot (WN/100, Q(1,:))
% set(gca, 'XDir', 'reverse')
% xlabel('Wavenumber [cm^-^1]')
% ylabel('Q (Scatter Efficiency)')
% title('Graph Showing Mie type scattering efficiency, Q')

end
