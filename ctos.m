function [LAT, LONG, RAD] = ctos(x,y,z)

%Conversion from Cartesian to Spherical Coordinates

LAT = atan (z/(sqrt(x^2+y^2)));
LONG = atan2 ( y , x );
RAD = sqrt ( x^2+y^2+z^2 );

if ~isnumeric (LAT) || ~isnumeric (LONG) || ~isnumeric (RAD)
    error ( 'Input arguments must be numeric. ');
end
if RAD <=0
    error ( 'Input argument r must be positive. ');
end
if any (size(LAT) ~= size(LONG)) || any (size(LAT) ~= size(RAD))
    error ('Latitude, Longtitude and Radius must be of the same size.');

end

