% Introduction to Satellite Geodesy - Exercise
% Assignment 3: Coordinate conversions and coordinate transformations

clc;
clear all; 
close all; 

% Exercise 1 --------------------------------------------------------------
%no code nedeed

% Exercise 2 --------------------------------------------------------------

% Inputes from Attachment 1:
  x=3782970.10; y=902154.92; z=5038375.59;
% Outputes:
% "ctos" - function for conversation - Cartesian to Spherical
% "stoc" - function for conversation - Spherical to Cartesian
[LAT,LONG,RAD]=ctos(x,y,z); % Cartesian to Spherical 1
[x1,y1,z1]=stoc(LAT,LONG,RAD); % Spherical to Cartesian 2
[LAT1,LONG1,RAD1]=ctos(x1,y1,z1); % Cartesian to Spherical 3

% Arclength of Spherical Coordinates;
[aLAT,aLONG]=arcle(LAT,LONG); [aLAT1,aLONG1]=arcle(LAT1,LONG1);

format long 
% Difference between Cartesian Coordinates 
C1=[x;y;z]; C2=[x1;y1;z1]; dC=C1-C2; %difference result

% Difference between Spherical Coordinates 
S1=[rad2deg(LAT); rad2deg(LONG); RAD]; 
S2=[rad2deg(LAT1); rad2deg(LONG1); RAD1]; format long 
dS=S1-S2; % difference result

% Difference between Spherical Coordinates (Arclength):
S3=[aLAT; aLONG]; S4=[aLAT1; aLONG1]; format long; dS1=S3-S4; 

% Exercise 3 --------------------------------------------------------------
% GRS80 Ellipsoid parameters: 
aG=6378137; fG=1/298.257222101; eG=sqrt(2*fG-fG^2); bG=aG*(1-fG);
% BESSEL Ellipsoid parameters: 
aB=6377397.155; fB=1/299.1528128; eB=sqrt(2*fB-fB^2); bB=aB*(1-fB);
% WGS84 Ellipsoid parameters: 
aW=6378137; fW=1/298.257223563; eW=sqrt(2*fW-fW^2); bW=aW*(1-fW);

% a) Conversation from Cartesian to Geodetic Coordinates
% GRS80 Ellipsoid:
[Bg,Lg,hg]=ctoge(x,y,z,fG,aG); % Cartesian to Geodetic
[xG1,yG1,zG1]=getoc(Bg,Lg,hg,fG,aG); % Geodetic to Cartesian  
[Bg1,Lg1,hg1]=ctoge(xG1,yG1,zG1,fG,aG); % Cartesian to Geodetic

% WGS84 Ellipsoid:
[Bw,Lw,hw]=ctoge(x,y,z,fW,aW); % Cartesian to Geodetic
[xW1,yW1,zW1]=getoc(Bw,Lw,hw,fW,aW); % Geodetic to Cartesian 
[Bw1,Lw1,hw1]=ctoge(xW1,yW1,zW1,fW,aW); % Cartesian to Geodetic

% Bessel Ellipsoid:
[Bb,Lb,hb]=ctoge(x,y,z,fB,aB); % Cartesian to Geodetic
[xB1,yB1,zB1]=getoc(Bb,Lb,hb,fB,aB); % Geodetic to Cartesian 
[Bb1,Lb1,hb1]=ctoge(xB1,yB1,zB1,fB,aB); % Cartesian to Geodetic

% Arc length with function distance
% GRS80
Bgd=rad2deg(Bg); arcBg=distance(0,0,Bgd,0,[aG,eG]);
Lgd=rad2deg(Lg); arcLg=distance(0,0,0,Lgd,[aG,eG]);
Bgd1=rad2deg(Bg1); arcBg1=distance(0,0,Bgd1,0,[aG,eG]);
Lgd1=rad2deg(Lg1); arcLg1=distance(0,0,0,Lgd1,[aG,eG]);
% WGS84W
Bwd=rad2deg(Bw); arcBw=distance(0,0,Bwd,0,[aW,eW]);
Lwd=rad2deg(Lw); arcLw=distance(0,0,0,Lwd,[aW,eW]);
Bwd1=rad2deg(Bw1); arcBw1=distance(0,0,Bwd1,0,[aW,eW]);
Lwd1=rad2deg(Lw1); arcLw1=distance(0,0,0,Lwd1,[aW,eW]);
% Bessel
Bbd=rad2deg(Bb); arcBb=distance(0,0,Bbd,0,[aB,eB]);
Lbd=rad2deg(Lb); arcLb=distance(0,0,0,Lbd,[aB,eB]);
Bbd1=rad2deg(Bb1); arcBb1=distance(0,0,Bbd1,0,[aB,eB]);
Lbd1=rad2deg(Lb1); arcLb1=distance(0,0,0,Lbd1,[aB,eB]);

% Arclength with approximation
% GRS80
sinBG=sin(Bg); e2G=2*fG-fG^2; NG=aG ./sqrt(1 - e2G * sinBG.^2);
BgA=Bg*sqrt(aG*bG); LgA=Lg*NG*cos(Bg);
sinBG1=sin(Bg1); NG1=aG ./sqrt(1 - e2G * sinBG1.^2);
Bg1A=Bg1*sqrt(aG*bG); Lg1A=Lg1*NG1*cos(Bg1);
% WGS84
sinBW=sin(Bw); e2W=2*fW-fW^2; NW=aW ./sqrt(1 - e2W * sinBW.^2);
BwA=Bw*sqrt(aW*bW); LwA=Lw*NW*cos(Bw);
sinBW1=sin(Bw1); NW1=aW ./sqrt(1 - e2W * sinBW1.^2);
Bw1A=Bw1*sqrt(aW*bW); Lw1A=Lw1*NW1*cos(Bw1);
% Bessel
sinBB=sin(Bb); e2B=2*fB-fB^2; NB=aB ./sqrt(1 - e2B * sinBB.^2);
BbA=Bb*sqrt(aB*bB); LbA=Lb*NB*cos(Bb);
sinBB1=sin(Bb1); NB1=aB ./sqrt(1 - e2B * sinBB1.^2);
Bb1A=Bb1*sqrt(aB*bB); Lb1A=Lb1*NB1*cos(Bb1);

% Difference between Cartesian Coordinates:
Cg1=[x;y;z]; Cg2=[xG1;yG1;zG1]; dCG=Cg1-Cg2; % GRS80
Cg3=[xW1;yW1;zW1]; dCG1=Cg1-Cg3; % WGS84
Cg4=[xB1;yB1;zB1]; dCG2=Cg1-Cg4; % Bessel

% Difference between Geodetic Coordinates (degree):
GE1=[rad2deg(Bg);rad2deg(Lg);hg];
GE2=[rad2deg(Bg1);rad2deg(Lg1);hg1]; dGE = GE1 - GE2; 
GE3=[rad2deg(Bw);rad2deg(Lw);hw];
GE4=[rad2deg(Bw1);rad2deg(Lw1);hw1]; dGE1 = GE3 - GE4; % WGS84
GE5=[rad2deg(Bb);rad2deg(Lb);hb];
GE6=[rad2deg(Bb1);rad2deg(Lb1);hb1]; dGE2 = GE5 - GE6; % Bessel

% Difference in arc length - distance function
GE8=[arcBg;arcLg]; GE9=[arcBg1;arcLg1]; dGE3=GE8-GE9; % GRS80
GE10=[arcBw;arcLw]; GE11=[arcBw1;arcLw1]; dGE4=GE10-GE11; % WGS84
GE12=[arcBb;arcLb]; GE13=[arcBb1;arcLb1]; dGE5=GE11-GE12; % Bessel

% Difference in arc length - simple approximation
GE14=[BgA;LgA]; GE15=[Bg1A;Lg1A]; dGE6=GE14-GE15; % GRS80
GE16=[BwA;LwA]; GE17=[Bw1A;Lw1A]; dGE7=GE16-GE17; % WGS84
GE18=[BbA;LbA]; GE19=[Bb1A;Lb1A]; dGE8=GE18-GE19; % Bessel

% Comparison in units ( degree, arclength )
COMP1 = GE4 - GE1; % in degree WGS84 - GRS80
COMP2 = GE10 - GE8; % in arc length WGS84 - GRS80 // distance function
COMP3 = GE16 - GE14; % in arc length WGS84 - GRS80 // simple aprrox.
COMP4 = GE5 - GE1; % in degree Bessel - GRS80
COMP5 = GE12 - GE8; % in arc length Bessel - GRS80 // distance function
COMP6 = GE18 - GE14; % in arc length Bessel - GRS80 //  simple approx.

% Comparison in percentage
COMP7 = (COMP1/GE1)*100; % in degree WGS84 - GRS80
COMP8 = (COMP2/GE8)*100; % arclength WGS84 - GRS80 // distance function
COMP9 = (COMP3/GE14)*100; % in arc length WGS84 - GRS80 // simple aprrox. 
COM10 = (COMP4/GE1)*100; % in degree Bessel - GRS80
COM11 = (COMP5/GE8)*100; % in arclength Bessel - GRS80 / distance function
COM12 = (COMP6/GE14)*100; % in arc length Bessel - GRS80 / simple approx.

% Ellipsoid globally closest to GRS80
GRS80 = [ aG; fG ] ;  WGS84 = [ aW; fW ]; Bessel = [ aB; fB ];
Gdiff = GRS80 - WGS84; Gdiff2 = (Gdiff).^2;
Gdiff1 = GRS80 - Bessel; Gdiff3 = (Gdiff1).^2;

if any ( Gdiff2 < Gdiff3 )
   disp ('WGS84 is globally closest to the GRS80.');
else 
   disp ('Bessel is globally closest to the GRS80.');
end
    
% Ellipsoid locally closest to GRS80
GRS80E=[Bg;Lg;hg]; WGS84E=[Bw;Lw;hw]; BESSE=[Bb;Lb;hb];
Ldiff = GRS80E - WGS84E; Ldiff2 = Ldiff .^2;
Ldiff1 = GRS80E - BESSE; Ldiff3 = Ldiff1 .^2;

if any ( Ldiff2 < Ldiff3 )
   disp ('WGS84 is locally closest to the GRS80.');
else 
   disp ('Bessel is locally closest to the GRS80.');
end

% Exercise 4 --------------------------------------------------------------
SPHall=[rad2deg(LAT);rad2deg(LONG);(RAD-6371); aLAT; aLONG ]; % SPHERE

GE80=[rad2deg(Bg);rad2deg(Lg);hg;arcBg;arcLg]; % GRS80
GE84=[rad2deg(Bw);rad2deg(Lw);hw;arcBw;arcLw]; % WGS84
GEBE=[rad2deg(Bb);rad2deg(Lb);hb;arcBb;arcLb]; % BESSEL

SEC1 = GE80 - SPHall; % GRS80 - SPHERE
SEcomp1 = sqrt( (SEC1(4,1)^2) + (SEC1(5,1)^2) + (SEC1(3,1)^2));
SEC2 = GE84 - SPHall; % WGS84 - SPHERE
SEcomp2 = sqrt( (SEC2(4,1)^2) + (SEC2(5,1)^2) + (SEC2(3,1)^2));
SEC3 = GEBE - SPHall; % BESSEL - SPHERE
SEcomp3 = sqrt( (SEC3(4,1)^2) + (SEC3(5,1)^2) + (SEC3(3,1)^2));

if ( SEcomp1 < SEcomp2 ) && ( SEcomp1 < SEcomp3 )
   disp ('GRS80 is locally closest to the sphere.');
elseif ( SEcomp2 < SEcomp1 ) &&  ( SEcomp2 < SEcomp3 )
   disp ('WGS84 is locally closest to the sphere.');
else 
   disp ('Bessel is locally closest to the sphere.');
end


% Exercise 5 --------------------------------------------------------------

% ITRF2008
x81=1492404.605; y81=(-4457266.520); z81=4296881.795; % 7205
x82=4075539.757; y82=931735.399; z82=4801629.449; % 7224
x83=5085442.772 ; y83=(2668263.635); z83=(-2768696.876); % 7232

vx81=(-0.0155); vy81=(-0.0012); vz81=0.0041; % 7205 velocities
vx82=(-0.0160); vy82=(0.0171); vz82=0.0101; % 7224 velocities
vx83=(-0.0015); vy83=(0.0196); vz83=0.0165; % 7232 velocities

% ITRF2005
x54=1492404.683; y54=(-4457266.515); z54=4296881.775; % 7205 
x55=4075539.836; y55=931735.313; z55=4801629.400; % 7224
x56=5085442.779; y56=2668263.544; z56=(-2768696.963); % 7232

% ep = 2005.0; % epoch of ITRF2008
% ep1 = 2000.0; % epoch of ITRF2005

% xme = ep + vx1 * ( ep - ep1 )
% yme = ep + vy1 * ( ep - ep1 ); 
% zme = ep + vz1 * ( ep - ep1 ); 

deltaX1 = [ (x54-x81); (y54-y81); (z54-z81) ];
deltaX2 = [ (x55-x82); (y55-y82); (y55-y82) ];
deltaX3 = [ (x56-x83); (y56-y83); (y56-y83) ];

l = [ deltaX1; deltaX2; deltaX3 ];

A = [ 1, 0, 0, x81, 0, z81, -y81;
      0, 1, 0, y81, -z81, 0, x81;
      0, 0, 1, z81, y81, -x81, 0; 
      1, 0, 0, x82, 0, z82, -y82;
      0, 1, 0, y82, -z82, 0, x82;
      0, 0, 1, z82, y82, -x82, 0; 
      1, 0, 0, x83, 0, z83, -y83;
      0, 1, 0, y83, -z83, 0, x83;
      0, 0, 1, z83, y83, -x83, 0; ]; 

%  P =  [ 1, 0, 0, 0, 0, 0, 0, 0, 0; 
%         0, 1, 0, 0, 0, 0, 0, 0, 0; 
%         0, 0, 1, 0, 0, 0, 0, 0, 0; 
%         0, 0, 0, 1, 0, 0, 0, 0, 0; 
%         0, 0, 0, 0, 1, 0, 0, 0, 0; 
%         0, 0, 0, 0, 0, 1, 0, 0, 0; 
%         0, 0, 0, 0, 0, 0, 1, 0, 0;  
%         0, 0, 0, 0, 0, 0, 0, 1, 0; 
%         0, 0, 0, 0, 0, 0, 0, 0, 1; ];  % Weight matrix P 9x9


% GS = (inv( A' * P * A )) * A' * P * l 
GS = (inv( A' * A )) * A' * l ; % without P

T = [ (GS(1,1)); (GS(2,1)); (GS(3,1))];
D = GS(4,1);
R = [ 0, - (GS(7,1)), (GS(6,1));
     (GS(7,1)), 0 , - (GS(5,1));
    -(GS(6,1)), (GS(5,1)), 0 ];

x87=-3950237.046; y87=2522347.621; z87=-4311562.205; % ITRF2008 7205
x58=-3950236.859; y58=2522347.586; z58=-4311562.417; % ITRF2005 7205 

Xt = x87 + T + D + R; 
Yt = y87 + T + D + R;   
Zt = z87 + T + D + R; 

xTransf = [ Xt(1,1); Yt(1,1); Zt(1,1)];
xGiven = [ x58; y58; z58 ]; 
xDIF = xTransf - xGiven;


