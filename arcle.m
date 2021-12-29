function [aLAT,aLONG] = arcle(LAT,LONG)
% Arclength of Spherical Coordinates
RE = 6371; %km     Earth Radius
aLAT = LAT * RE; % LONG = constant
aLONG = LONG * RE * cos (LAT); % LAT = constant
end

