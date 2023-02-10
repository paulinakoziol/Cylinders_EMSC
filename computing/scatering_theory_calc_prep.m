function [ZDis] = scatering_theory_calc_prep(ZRef, r_min, r_max, n_min, n_max, ...
    spacings, zeta, orientation)
% This function performs Kramers Kronig transformation and prepares
% permutations of parameters (radius and refractive index) for further
% Scattering efficiency calculations 

% More info about naming convension and physical background availible in:
% Analyst, 2010, 135, 268–277

WN = str2num(ZRef.v)';

%% Kramers-Kronig transform
ref_n = kkre(WN, ZRef.d);
ref_n = ref_n/abs(min(ref_n));

%% Combination of parameters r, a, b
rr = 1e-6*linspace(r_min, r_max, spacings); %radius range with conversion from um to m
aa = linspace(n_min, n_max, spacings); %refractive index range
bs = spacings;

number_rows = length(rr)*length(aa)*bs;
% Permutations of parameters:
z = zeros(number_rows, 3);
j = 0;
for k = 1 : length(rr)
    r = rr(k);
    for n = 1 : length(aa)
        a = aa(n);
        bb = linspace(0.0, a-1.01, bs);
        for i = 1 : length(bb)
            j = j + 1;
            b = bb(i);
            z(j,1:3) = [r, a, b];
        end
    end
end

% Real part of the refractive index as described by equation 11 of 
% Analyst, 2010, 135, 268–277
ri = z(:,2) + z(:,3)*ref_n; 

%% Main function for scattering efficiency calculation
out = full_theory_cylinders(WN, z(:,1), ri, zeta, orientation);

ZDis.d = out;
ZDis.v = ZRef.v;
ZDis.i = num2str ([1:number_rows]');

end






