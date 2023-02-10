function [Q] = scattering_efficiency(x, A, B, C, D, V, W, orientation)
% scattering_efficiency function calculates scattering efficiency for
% samples with cylindrical domains using full Mie type theory provided in:
% Applied Spectroscopy 2019, Vol. 73(8) 859–869

% Denominator calculation
denom = (W.*V+1i*(D).^2);

% a_l and b_l coefficients calculation
a_l_lon = (C.*V-B.*D)./denom;
a_l_tra = -(A.*V-1i*C.*D)./denom;
b_l_lon = (W.*B+1i*D.*C)./denom;
b_l_tra = -1i*(C.*W+A.*D)./denom;

% Final scattering efficiency calculation considering polarization
if orientation == -1 
    % Unpolarized light case
    Q_lon = (2./x).*(abs(b_l_lon(:,:,1)).^2+2*sum((abs(b_l_lon(:,:,2:end)) ...
        .^2+abs(a_l_lon(:,:,2:end)).^2), 3));
    Q_tra = (2./x).*(abs(a_l_tra(:,:,1)).^2+2*sum((abs(a_l_tra(:,:,2:end)) ...
        .^2+abs(b_l_tra(:,:,2:end)).^2), 3));
    Q = 0.5*(Q_lon+Q_tra);

elseif orientation == 0  
    % Polarization paralell to the cylinder axis
    Q = (2./x).*(abs(b_l_lon(:,:,1)).^2+2*sum((abs(b_l_lon(:,:,2:end)) ...
        .^2+abs(a_l_lon(:,:,2:end)).^2), 3));
    
elseif orientation == 90  
    % Polarization perpendicular to the cylinder axis
    Q = (2./x).*(abs(a_l_tra(:,:,1)).^2+2*sum((abs(a_l_tra(:,:,2:end)) ...
        .^2+abs(b_l_tra(:,:,2:end)).^2), 3));

else  
    orientation = orientation*pi/180;
    % Other linear polarization cases
    Q_lon = (2./x).*(abs(b_l_lon(:,:,1)).^2+2*sum((abs(b_l_lon(:,:,2:end)) ...
        .^2+abs(a_l_lon(:,:,2:end)).^2), 3));
    Q_tra = (2./x).*(abs(a_l_tra(:,:,1)).^2+2*sum((abs(a_l_tra(:,:,2:end)) ...
        .^2+abs(b_l_tra(:,:,2:end)).^2), 3));
    Q=Q_lon*cos(orientation)^2+Q_tra*sin(orientation)^2;
end
    
end

