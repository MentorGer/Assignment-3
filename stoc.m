function [x, y, z] = stoc( LAT, LONG, RAD)

% Conversation from Spherical Coordinates to Cartesian 

x = RAD * cos (LONG)* cos (LAT);
y = RAD * sin (LONG)* cos (LAT);
z = RAD * sin(LAT);

end

