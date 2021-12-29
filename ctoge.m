function [B,L,h] = ctoge(x,y,z,f,a)
% Method of BOWRING
% Conversion from cartesian coordinates to geodetic coordinates
e2 = 2 * f - f ^ 2;
ep2 = e2 / (1 - e2); % square of second eccentricity
b = a * (1 - f);
L = atan2(y,x); % geodetic longitude
p = sqrt(abs(x).^2+abs(y).^2);% distance from z-axis
beta = atan2(z, (1 - f) * p); % parametric latitude start value
B = atan2(z + b * ep2 * sin(beta).^3,...
p - a * e2 * cos(beta).^3);
betaNew = atan2((1 - f)*sin(B), cos(B));
i_ter = 0;
while any(beta(:) ~= betaNew(:)) && i_ter < 5
    beta = betaNew;
    B = atan2(z + b * ep2 * sin(beta).^3,...
    p - a * e2 * cos(beta).^3);
    betaNew = atan2((1 - f)*sin(B), cos(B));
    i_ter = i_ter + 1;
end
% calculate ellipsoidal height from the final value of latitude
sinB = sin(B);
N = a ./ sqrt(1 - e2 * sinB.^2);
h = p .* cos(B) + (z + e2 * N .* sinB) .* sinB - N;
end

