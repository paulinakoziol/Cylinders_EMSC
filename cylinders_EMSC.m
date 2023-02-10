function [WN_out, ZCorr, history, timing] = cylinders_EMSC(WN, ZRaw, ...
    correction_options, WN_Ref_in, Ref_in)

% cylinders_EMSC Resonant Mie type Scattering Correction for samples with
% cylindrical domains
% 
% Syntax
%   [WN_out, ZCorr, history] =
%       cylinders_EMSC(WN, ZRaw, correction_options);
%   [WN_out, ZCorr, history] =
%       cylinders_EMSC(WN, ZRaw, correction_options, WN_Ref_in, Ref_in);
%
% Description
%   [WN_out, ZCorr, history] = cylinders_EMSC(WN, ZRaw, correction_options,orientation)
%   performs a correction of mid-infrared spectra corresponding to resonant
%   Mie type theory for cylindrical samples. WN is the wavenumber vector. ZRaw
%   is a 2D matrix containing absorbance spectra in rows. correction_options 
%   is a vector containing options (see below). WN_out is the wavenumber 
%   vector of the corrected spectra. ZCorr is a 2D matrix of corrected spectra in rows. history is
%   a 3D matrix of dimension; numberOfSpectra x length(WN_out) x
%   iterations. timing is a vector with calcuation time for each spectrum.
%
%   [WN_out,ZCorr,history] =
%   cylinders_EMSC(WN,ZRaw,correction_options,WN_Ref_in,Ref_in) utilises a
%   user supplied reference spectrum. WN_Ref_in is a wavenumber vector
%   corresponding to the reference spectrum. Ref_in is a vector of absorbance
%   values of the reference sample.
%
%   correction options is a vector with the following default parameters:
%   correction_options = [ ...
%     1000 ; % 1. Lower wavenumber range
%     4000 ; % 2. Upper wavenumber range
%     10   ; % 3. Number of iterations
%     7    ; % 4. Number of principal components used (7 recommended)
%     2    ; % 5. Lower range for scattering particle diameter / um
%     8    ; % 6. Upper range for scattering particle diamter / um
%     1.1  ; % 7. Lower range for average refractive index
%     1.5  ; % 8. Upper range for average refractive index
%     10   ; % 9. Number of values for each scattering parameter default 10
%     1    ; % 10. Orthogonalisation, 0 = no, 1 = yes. (1 recommended)
%     1    ; % 11. Which reference spectrum, 1 = Matrigel, 2 = PCL spectrum , 3 = Mean spectrum
%     pi/2 ; % 12. Incidence angle (zeta) (perpendicular - 90deg = pi/2)
%     -1  ]; % 13. Cylinder orientation with respect to linear polarization
%
%   The correction_options must be supplied by the user, defaults are not
%   provided by this function. 


%% Correction options

lower_WN = correction_options(1);
upper_WN = correction_options(2);
iterations = correction_options(3);
NCOMP = correction_options(4);
r_min = correction_options(5);
r_max = correction_options(6);
n_min = correction_options(7);
n_max = correction_options(8);
spacings = correction_options(9);
GSP_flag = correction_options(10);
ref_option = correction_options(11);
zeta = correction_options(12);
orientation = correction_options(13);

%% Wavenumber vector adjustment

WN = make_column(WN);

[a, b] = find_value_min_max(WN, lower_WN, upper_WN);

WN = WN(a:b);
ZRaw = ZRaw(:,a:b);

[N, K] = size(ZRaw);
spectra_names = num2str((1:N)');

%% Reference spectrum adjustment

if (nargin == 3)
    if (ref_option == 1)
        load Matrigel_Reference_Raw
        ZRef = spline(ZRef_Raw(:,1), ZRef_Raw(:,2), WN)'; 
        clear ZRef_Raw
    elseif (ref_option == 2)
        load PCL_ZRef_Raw; 
        ZRef = spline(ZRef_Raw(:,1), ZRef_Raw(:,2), WN)';
        clear ZRef_Raw
    elseif (ref_option == 3)
        ZRef = mean(ZRaw,1);
    else
        error('No reference spectrum option selected')
    end
end

if (nargin == 5)
    if ((length(WN) == length(WN_Ref_in)) ...
            && (length(ZRaw) == length(Ref_in)))
        ref_length_check = 1;
    else
        ref_length_check = 0;
    end
    
    if (ref_length_check == 0)
        ZRef = spline(WN_Ref_in, Ref_in, WN)';
    else
        ZRef = Ref_in;
    end
end

ZRef = ZRef/max(ZRef); %Reference spectrum normalization

%% Make Saisir Structures & Weights

temp = ZRaw; 
clear ZRaw;

ZRaw.d = temp; 
clear temp;

ZRaw.v = num2str(WN);
ZRaw.i = spectra_names;

temp = ZRef; 
clear ZRef

ZRef.d = temp;
ZRef.v = ZRaw.v;
ZRef.i = 'Ref';

ZWeightsSpec = Spectrum_Weights(ZRaw);

if (iterations > 1)
    history = zeros(N, K, iterations);    
end

ZWeightsSpec = Down_weight_spectrum(ZWeightsSpec, 2300, 2400, 0.01); %CO2 weights go down

%% Correction using cylinders-EMSC for samples with cylindrical domains

startT=tic;
[ZCorr, mod_para] = scattering_maths(ZRef, ZRaw, ZWeightsSpec, ...
    NCOMP, GSP_flag, r_min, r_max, n_min, n_max, spacings, zeta, orientation);
history(:,:,1) = ZCorr.d;
res_hist(:,:,1) = mod_para.residual; %#ok<NASGU>

timing=zeros(1,N);

if (iterations > 1)
    for j = 1:N
        ZRaw2.d = ZRaw.d(j,:);
        for i = 2:iterations
            try %#ok<TRYNC>
                ZRef.d = history(j,:,i-1);
                ZRef.d = fit_gauss_whole_spec_split(WN, ZRef.d')';
                
                [ZCorr, mod_para] = scattering_maths(ZRef, ZRaw2, ZWeightsSpec, ...
                    NCOMP, GSP_flag, r_min, r_max, n_min, n_max, ...
                    spacings, zeta, orientation); %#ok<ASGLU>
                
                history(j,:,i) = ZCorr.d;

%                 disp(['Spectrum ', num2str(j), ' Iteration ', num2str(i), '  ', datestr(now)])
            end
        end
        ZCorr.d = history(:,:,end);
        timing(1,j) = toc(startT);
    end
end

ZCorr.d = history(:,:,end);

WN_out = WN;
ZCorr = ZCorr.d;

end % end of function cylinders_EMSC
