function [q, mod_par] = scattering_maths(R, S, W, NCOMP, GSP_flag, r_min, ...
    r_max, n_min, n_max, spacings, zeta, orientation)

q = S;
[N K] = size(S);

%% Scattering efficiency for a current iteration
[ZDis] = scatering_theory_calc_prep(R, r_min, r_max, n_min, n_max, spacings, ...
    zeta, orientation);

%% PCA decomposition and model building
[T, P] = NIPALS_easy(ZDis.d, NCOMP); 

model(:,1) = R.d';
model(:,2) = ones(length(R.d), 1);
model(:,3) = linspace(0, 1, length(R.d));
model(:,4:3+NCOMP) = P(:,1:NCOMP);

%% Gram-Schimd Process
if GSP_flag == 1
    model = GSP(model, 101);
end

%% Least Squares Part
cons = lscov(model, S.d', W.d');

a = cons(2:end,:);
b = model(:,2:end)';

scatter = a'*b;

mod_par.cons = cons;
mod_par.model = model;
mod_par.estimated = cons'*model';
mod_par.residual = S.d - mod_par.estimated;
mod_par.scatter = scatter;

q.d = S.d - scatter;

end % end of Mie_maths



























