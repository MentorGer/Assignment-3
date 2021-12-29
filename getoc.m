function [x,y,z] = getoc(B,L,h,f,a)
%Conversation from geodetic coordinates to cartesian coordinates
e2 = 2 * f - f ^ 2;
ep2 = e2 / (1 - e2); % square of second eccentricity
b = a * (1 - f);
sinB = sin(B);
N = a ./ sqrt(1 - e2 * sinB.^2);
x = (N+h)*cos(B)*cos(L);
y = (N+h)*cos(B)*sin(L);
z = (N-e2*N+h)*sin(B);
end

