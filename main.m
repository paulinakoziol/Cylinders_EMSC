clear all
close all

% This code requirec GPU unit

addpath(genpath('computing'));,
addpath(genpath('Model_data'));

file = 'PCL_fiber.mat'; % Data for thin PCL fiber measured with unpolarized light 

%% Load and arrange needed data 
load('mask.mat');
load('frequencies.mat');
load(file)

[A, B, numVar] = size(data);
figure; imshow(data(:,:,237));

%Excluding background pixels
dataR = reshape(data, A*B, numVar);
mask = reshape(mask, A*B, 1);
dataR = dataR(mask,:); 
figure; plot(frequencies,dataR(2,:));

%% Correction options
    
correction_options = [ ...
    898  ;      % 1. Lower wavenumber range
    3845 ;      % 2. Upper wavenumber range
    10   ;      % 3. Number of iterations
    7    ;      % 4. Number of principal components used (7 default)
    5    ;      % 5. Lower range for scattering particle radius / um
    15   ;     % 6. Upper range for scattering particle radius / um
    1.3  ;      % 7. Lower range for average refractive index
    1.7  ;      % 8. Upper range for average refractive index
    10   ;      % 9. Number of values for each scattering parameter (a,b,d) default 10
    1    ;      % 10. Gram-Schmidt Process option, 0 = no, 1 = yes. (1 recommended)
    2    ;      % 11. Which reference spectrum, 1 = Matrigel, 2 = PCL spectrum , 3 = Mean spectrum
    pi/2 ;      % 12. Incidence angle (zeta) (perpendicular - 90deg = pi/2)
    -1  ];      % 13. Cylinder orientation with respect to linear polarization (degrees):
                %     from 0 for paralell to 90 for perpendicular orientation,
                %     -1 for unpolarized light.                   

%% Data preparation
dataR = fliplr(dataR); % Spectra to be corrected, one spectrum per row
WN = flipud(frequencies); % Wavenumbers, coloumn vector 
    
N = size(dataR, 1);  %Number of spectra for correction
disp(['There is ' num2str(N) ' spectra for correction...'])
    
%% The actual correction
[WN_corr, correctedSpectra, history, timing] = cylinders_EMSC(WN, ...
    dataR(1:3,:), correction_options); %Examplary correction for three pixels
   
figure; plot(WN_corr,correctedSpectra(1,:)); hold on; plot(WN, dataR(1,:))

% %% Shaping corrected data back into the 3D structure 
% % only possible if the full dataset is corrected
% 
% correctedSpectra = fliplr(correctedSpectra);
% 
% dataC = zeros(A*B, numVar);
% dataC(mask,:) = correctedSpectra; 
% dataC = reshape(dataC, A, B,numVar);
% figure; imshow(dataC(:,:,237));
% 


