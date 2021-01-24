% *********************************************************************
% Aperture.m
% *********************************************************************
% Programmed by Xin Xie, 11/2002
% Reference: "Antenna Theory: Analysis and Design", Third Edition '
% Constantine A. Balanis, Chapter 13, Aperture Antennas
% ---------------------------------------------------------------------
% The program calculates the radiation characteristics of an aperture
% antenna either:
% a. Mounted on an infinite perfect electric conducting (PEC) ground plane
% b. Not mounted on a PEC ground plane
%
% The aperture antennas considered are either:
% a. Rectangular of dimensions a and b
% b. Circular of radius a
%
% Each aperture has either of the following aperture field distributions:
% a. Uniform
% b. Dominant mode:
%    1.  TE10 for the rectangular
%    2.  TE11 for the circular
%------------------------------------------------------------------------
% The radiation characteristics calculated are based on theformulations 
% of Chapter 12, Sections 12.3-12.6, and Tables 12.1 and 12.2. 
%-------------------------------------------------------------------------
% For each aperture, the following radiation characteristics are specified:
%
% INPUT PARAMETERS:
%
% Aperture Antenna Type:
% a. Rectangular
% b. Circular
%
% Aperture Dimensions
% a. Rectangular
%    1. a (in lambda)
%    2. b (in lambda)
% b. Circular
%    1. Radius a (in lambda)
%
% Aperture field distribution:
% a. Rectangular
%    1. Uniform
%    2. Dominant Mode (TE10)
% b. Circular
%    1. Uniform
%    2. Dominant Mode (TE11)
% Aperture mounting
% a. Infinite Perfect Electric Conductor (PEC)
% b. No Ground Plane
%----------------------------------------------------------------
%OUTPUT PARAMETERS:
%----------------------------------------------------------------
% Aperture efficiency
%---------------------------------------------------------------
% Approximate Directivity D0 - Table 12.1
% Approximate Directivity D0 (dB) - Table 12.1
% Approximate HPBW  for E-plane (in degrees) - Table 12.1
% Approximate HPBW  for H-plane (in degrees) - Table 12.1 
% Approximate FNBW  for E-plane (in degrees) - Table 12.1 
% Approximate FNBW  for H-plane (in degrees) - Table 12.1 
%----------------------------------------------------------------
% Radiated power (watts) - Integration
% Numerical Directivity D0 (dimensionless) - Integration
% Numerical HPBW  for E-plane (in degrees) - Integration
% Numerical HPBW  for H-plane (in degrees) - Integration
% Numerical FNBW  for E-plane (in degrees) - Integration
% Numerical FNBW  for H-plane (in degrees) - Integration
%----------------------------------------------------------------
% Plots for both rectangular and circular apertures:
% a. E-Plane normalized far-field amplitude pattern (in dB)
% b. H-plane normalized far-field amplitude pattern (in dB)
% c. E-plane directivity pattern (in dB)
% d. H-plane directivity pattern (in dB)
%----------------------------------------------------------------
% The program uses 'polar_dB" to plot the polar patterns.
% ---------------------------------------------------------------
% To get the correct output, follow the instructions below: 
%----------------------------------------------------------------
%	Input
%		Choice:
%			A: Rectangular (Dimensions (a,b) (in lamda))
%			B: Circular (a (in lamda))
%		Choice: Distribution
%			A: Uniform
%			B: Dominant mode(Rect: TE10; Circular: TE11)
%		Choice:
%			A: Ground plane
%			B: No ground plane

%	Output
%		Directivity (dimensionless in dB)
%			1. Approximate
%			2. Numerical (directivity program)
%		Aperture efficiency (%)
%		HPBW (in degrees)
%			1. E-plane
%				a. approximate
%				b. numerical (iterative)
%			2. H-plane
%				a. approximate
%				b. numerical (iterative)
%		FNBW (in degrees)
%			1. E-plane
%				a. approximate
%				b. numerical (iterative)
%			2. H-plane
%				a. approximate
%				b. numerical (iterative)
%		PLOT: (Polar: in dB)
%			1. Normalized amplitude pattern (0-60 dB)
%			2. Directivity pattern (60 dB dynamic range)
function aperture

clear all;
close all;

format long;
global atype adis agr aca ara arb;

rtod = 180/pi;
err = 0;
zeroj01 = 2.4048;
zeroj11 = 3.8317;


%---Choice of output---

fprintf('Output device option \n\tOption (1):  Screen\n\tOption (2):  File \n');
ERR = 1;
while(ERR ~= 0)
   DEVICE = str2num(input('Output device = ','s'));
   if(DEVICE == 1)
      ERR = 0;
   elseif(DEVICE == 2)
      FILNAM = input('Input the desired output filename: ','s');
      ERR = 0;
   else
      fprintf('\nOutputting device number should be either 1 or 2\n');
   end
end



fprintf('\n\nPlease choose Aperture Antenna Type:\n\t(1).  Rectangular\n\t(2).  Circular\n');
atype = input('Enter the Option (1 or 2): ');

strtype = 'Aperture Antenna: ';

if (atype == 1)
   fprintf('\n\nPlease enter aperture dimensions a,b (in lambda)\n');
   fprintf('\t');
   ara = input('a (in lambda) = ');
   fprintf('\t');
   arb = input('b (in lambda) = ');
   
   strtype = 'Rectangular';
   strdimen = strcat('a (in lambda)= ',num2str(ara),', b (in lambda)= ',num2str(arb));
      
elseif (atype == 2)
   fprintf('\n\nPlease enter the radius a (in lambda)\n\t');
   aca = input('a (in lambda) = ');
   
   strtype = 'Circular';
   strdimen = strcat('a (in lambda)= ',num2str(aca));
   
else
   fprintf(' The input is wrong!!!\n');
   err = 1;
   return;
end

fprintf('\n\nPlease choose the E-field distribution: \n\t(1).  Uniform\n\t(2).  Dominant mode( Rect: TE10; Circ: TE11)\n');
adis = input('Enter the Option (1 or 2): ');
if adis ~= 1 & adis ~= 2
   err = 1;
end

if adis == 1
   strdis = 'Uniform';
else
   if atype == 1
      strdis = 'TE10';
   else
      strdis = 'TE11';
   end
end

if err == 1
   fprintf(' The input is wrong!!!\n');
   return;
end

fprintf('\n\nPlease choose "Ground plane" option: \n\t(1).  Ground plane\n\t(2).  No ground plane\n');
agr = input('Enter the Option (1 or 2): ');
if agr ~= 1 & agr ~= 2
   err = 1;
end

if err == 1
   fprintf(' The input is wrong!!!\n');
	return;   
end

if agr == 1
   strgrd = 'With ground plane';
else
   strgrd = 'No ground plane';
end

% if draw 3D pattern?
d3d = 0;
%fprintf('\n');
%strd3d = input('If the 3D amplitude patterns are drawn or not? (yes/no): ', 's');
%if strd3d == 'y'
%   d3d = 1;
%end

% avoid sin(0)/0
shift = 0.01;

% rectangular aperture
if (atype == 1)
   % aperture efficiency
   aeff = 1.0;
   
   % approximate HPBW
   phpbwE = 50.6/arb;
   phpbwH = 50.6/ara;
   % approximate FNBW
   pfnbwE = 114.6/arb;
   pfnbwH = 114.6/ara;
   
   if (adis ==1 & agr ==1)
      % Approximate directivity
      pD0 = 4*pi*ara*arb;
      
      % numerical directivity
      nfnbwE = asin(1/arb)*rtod*2;
      nfnbwH = asin(1/ara)*rtod*2;
      %nhpbwE = asin(2.782/(2*pi*arb))*rtod*2;
      %nhpbwH = asin(2.782/(2*pi*ara))*rtod*2;

      
         
   elseif(adis ==1 & agr==2)
      % Approximate directivity
      pD0 = 4*pi*ara*arb;
      
      % numerical
      nfnbwE = asin(1/arb)*rtod*2;
      nfnbwH = asin(1/ara)*rtod*2;

      
   elseif(adis == 2 & agr == 1)
      % Approximate directivity
      pD0 = 0.81*4*pi*ara*arb;
      aeff = 0.81;
      
      % Approximate HPBW and FNBW in H-plane
      phpbwH = 68.8/ara;
      pfnbwH = 171.9/ara;
      
      % numerical 
      nfnbwE = asin(1/arb)*rtod*2;
      % nfnbwH = asin(1/(2*ara))*rtod*2;
      rtheta = asin(2/ara);
		% fprintf('\nComputing the FNBW in H-plane ... ');
      nfnbwH = bsect(3,shift,rtheta-shift);
      nfnbwH = nfnbwH*rtod*2;
      
   elseif(adis == 2 & agr == 2)
      % Approximate directivity
      pD0 = 0.81*4*pi*ara*arb;
      aeff = 0.81;
      
      % Approximate HPBW and FNBW in H-plane
      phpbwH = 68.8/ara;
      pfnbwH = 171.9/ara;
      
      % numerical 
      nfnbwE = asin(1/arb)*rtod*2;
      % nfnbwH = asin(1/(2*ara))*rtod*2;
      rtheta = asin(2/ara);
		% fprintf('\nComputing the FNBW in H-plane ... ');
      nfnbwH = bsect(3,shift,rtheta-shift);
      nfnbwH = nfnbwH*rtod*2;

   else
      fprintf('This situation can''t be computed by this program!!!');
      return;
   end
   
   % fprintf('\nComputing the HPBW in E-plane ... ');
   nhpbwE = bsect(1,shift,nfnbwE/(2*rtod)-shift);
   nhpbwE = nhpbwE*rtod*2;
   % fprintf('\nComputing the HPBW in H-plane ... ');
   nhpbwH = bsect(2,shift,nfnbwH/(2*rtod)-shift);
   nhpbwH = nhpbwH*rtod*2;
   
   pD0dB = 10*log10(pD0);
      
   
% circular aperture   
else
   aeff = 1.0;
   phpbwE = 29.2/aca;
   pfnbwE = 69.9/aca;
   
   if (adis == 1) % & agr == 1)
      pD0 = (2*pi*aca)^2;
      phpbwH = phpbwE;
      pfnbwH = pfnbwE;
      
      % numerical
      nfnbwE = asin(zeroj11/(2*pi*aca))*rtod*2;
      nfnbwH = asin(zeroj11/(2*pi*aca))*rtod*2;
      
   elseif (adis == 2) % & agr == 1)
      aeff = 0.836;
      pD0 = aeff*(2*pi*aca)^2;

      phpbwH = 37.0/aca;
      pfnbwH = 98.0/aca;
      
      % numerical
      nfnbwE = asin(zeroj11/(2*pi*aca))*rtod*2;
      % fprintf('\nComputing the FNBW in H-plane ... ');
      
      rtheta = asin(1/aca);
      % rtheta = 1;
      nfnbwH = bsect(3,shift,rtheta);
      nfnbwH = nfnbwH*rtod*2;
      
   else
      fprintf('This situation can''t be computed by this program!!!');
      return;
   end
   
   % fprintf('\nComputing the HPBW in E-plane ... ');   
   nhpbwE = bsect(1,shift,nfnbwE/(2*rtod)-shift);
   nhpbwE = nhpbwE*rtod*2;
   
   % fprintf('\nComputing the HPBW in H-plane ... ');
   nhpbwH = bsect(2,shift,nfnbwH/(2*rtod)-shift);
   nhpbwH = nhpbwH*rtod*2;
   
   pD0dB = 10*log10(pD0);
   
end

fprintf('\n\nThe computation will take about one minute, please wait ... \n');

% numerical for the directvity (directivity program)
%---Definition of constants and initialization----

PI = 4.0*atan(1.0);
PRAD = 0.0;

UMAX = U(0,0);
UMIN = UMAX;


%---Read lower and upper bounds on theta and phi---

TL = 0;

if agr == 1
   TU = 90;
else
   TU = 180;
end

PL = 0;
PU = 360;

%---Main Program ---

THETA = PI/180.0;
PHI = PI/180.0;
PLL = PL+1;

while(PLL <= PU)
   XJ = PLL*PI/180.0;
   TLL = TL+1;
   while(TLL <= TU)
      XI = TLL*PI/180.0;
      F = U(XI,XJ);
      if(F > UMAX)
         UMAX = F;
      end
      if(F < UMIN)
         UMIN = F;
      end
      
      UA = THETA*PHI*F*sin(XI);
      PRAD = PRAD+UA;
      TLL = TLL+1;
   end
   PLL = PLL+1;
end

DR = 4.0*PI*UMAX/PRAD;
DRDB = 10.0*log10(DR);

% to plot the amplitude pattern
fprintf('\nPlotting the amplitude and directivity pattern ...\n');

stheta = 0.1;
rmin = -60;
rmax = 0;
rticks = 4;
linestyle ='-';

% ezsurf('cos(s)*cos(t)','cos(s)*sin(t)','X(s,t)',[0,pi/2,0,2*pi])
% E plane phi = pi/2
i = 1;

if agr == 1
   sdeg = -90;
   edeg = 90;
else
   sdeg = 0;
   edeg = 360;
end

for ii = sdeg:stheta:edeg
   theta(i) = ii/rtod;
   de(i) = amp(1,theta(i));
   dh(i) = amp(2,theta(i));
   if ii == 0
      denorm = de(i);
      dhnorm = dh(i);
   end
   
   i = i+1;
end

% normalize the amplitude pattern
de = de/denorm;
for i = 1:length(de)
   
   if de(i) == 0
      de(i) = rmin;
   else
      de(i) = 10*log10(de(i));
   end
   
   if de(i) < rmin
      de(i) = rmin;
   end
   
end

dh = dh/dhnorm;
for i = 1:length(dh)
   if dh(i) == 0
      dh(i) = rmin;
   else
      dh(i) = 10*log10(dh(i));
   end
   
   if dh(i) < rmin
      dh(i) = rmin;
   end
   
end

figure(1);
polar_dB(theta,de,rmin,rmax,rticks,linestyle);
title('Normalized amplitude pattern in \itE-plane\rm (in dB)');

figure(2);
polar_dB(theta,dh,rmin,rmax,rticks,linestyle);
title('Normalized amplitude pattern in \itH-plane\rm (in dB)');


% to plot directivity 3D pattern
%if d3d == 1
%   figure(5);
%   z = cplxgrid(40);
%   cplxmap(z,func3d(z)./denorm);
%end

% to plot the directivity pattern

% E-plane
% function hpol = polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%for theta = 0:1:360
i = 1;
for ii = sdeg:stheta:edeg
   
   dtheta(i) = ii/rtod;
   ded(i) = U(dtheta(i),pi/2);
   dhd(i) = U(dtheta(i),0);
   i = i+1;
end

ded = abs(ded);
ded = 4.0*pi*ded./PRAD;
for i = 1:length(ded)
   
   if ded(i) == 0
      ded(i) = rmin;
   else
      ded(i) = 10*log10(ded(i));
   end
   
   if ded(i) < rmin
      ded(i) = rmin;
   end
   
end

dhd = abs(dhd);
dhd = 4.0*pi*dhd./PRAD;
for i = 1:length(dhd)
   
   if dhd(i) == 0
      dhd(i) = rmin;
   else
      dhd(i) = 10*log10(dhd(i));
   end
   
   if dhd(i) < rmin
      dhd(i) = rmin;
   end
   
end


figure(3);
rmax = DRDB;
%rmin = 4.0*pi*UMIN/PRAD;
%rmin = 10.0*log10(rmin);

polar_dB(dtheta,ded,rmin,rmax,rticks,linestyle);
title('Directivity pattern in \itE-plane\rm (in dB)');

figure(4);
polar_dB(dtheta,dhd,rmin,rmax,rticks,linestyle);
title('Directivity pattern in \itH-plane\rm (in dB)');



%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'wt');
else
   fid = DEVICE;
end


fprintf(fid,'\n\n');

fprintf(fid,'APERTURE:  '); fprintf(fid,strtype);

fprintf(fid,'\n----------------------------------------------------------------');
fprintf(fid,'\nInput parameters:\n----------------------------------------------------------------\n');
fprintf(fid,'Aperture Antenna Type:  '); fprintf(fid, strtype);
fprintf(fid,'\nAperture Dimensions:  '); fprintf(fid, strdimen);
fprintf(fid,'\nField Distribution:  '); fprintf(fid, strdis);
fprintf(fid,'\n'); fprintf(fid, strgrd);

fprintf(fid, '\n----------------------------------------------------------------');

fprintf(fid, '\nOutput parameters:\n----------------------------------------------------------------');
fprintf(fid, '\nThe aperture efficiency is = %3.4f\n', aeff);

fprintf(fid, '\nApproximate Directivity D0 (dimensionless) = %3.4f ', pD0);
fprintf(fid, '\nApproximate Directivity D0 (dB)            = %3.4f ', pD0dB);
fprintf(fid, '\nApproximate HPBW  for E-plane (in degrees) = %3.4f ', phpbwE);
fprintf(fid, '\nApproximate HPBW  for H-plane (in degrees) = %3.4f ', phpbwH);
fprintf(fid, '\nApproximate FNBW  for E-plane (in degrees) = %3.4f ', pfnbwE);
fprintf(fid, '\nApproximate FNBW  for H-plane (in degrees) = %3.4f ', pfnbwH);

fprintf(fid, '\n----------------------------------------------------------------');

%---Echo input parameters and output computed parameters---
fprintf(fid, '\nRadiated power (watts) = %6.4f',PRAD);
fprintf(fid, '\nNumerical Directivity D0 (dimensionless)   = %3.4f',DR);
fprintf(fid, '\nNumerical Directivity D0 (dB)              = %3.4f',DRDB);
%---End of program-------------------------------------------

fprintf(fid, '\nNumerical HPBW  for E-plane (in degrees)   = %3.4f ', nhpbwE);
fprintf(fid, '\nNumerical HPBW  for H-plane (in degrees)   = %3.4f ', nhpbwH);
fprintf(fid, '\nNumerical FNBW  for E-plane (in degrees)   = %3.4f ', nfnbwE);
fprintf(fid, '\nNumerical FNBW  for H-plane (in degrees)   = %3.4f ', nfnbwH);



fprintf(fid,'\n');

fclose(fid);



function y = U(theta, phi)
global atype adis agr ara arb aca;

if atype == 1
   X = pi*ara*sin(theta)*cos(phi);
   if X == 0
      sincX = 1;
   else
      sincX = sin(X)/X;
   end
   
   Y = pi*arb*sin(theta)*sin(phi);
   if Y == 0
      sincY = 1;
   else
      sincY = sin(Y)/Y;
   end
   
   if (adis == 1 & agr == 1)
      y = (sin(phi)^2+(cos(theta)*cos(phi))^2)* ...
         (sincX*sincY)^2;
   elseif (adis == 1 & agr == 2)
      y = ((1+cos(theta))/2*sincX*sincY)^2;
   elseif (adis == 2 & agr == 1)
      y = (sin(phi)^2+(cos(theta)*cos(phi))^2)* ...
         (cos(X)/(X^2-(pi/2)^2)*sincY)^2;
   elseif (adis == 2 & agr == 2)
      y = ((1+cos(theta))/2*cos(X)/(X^2-(pi/2)^2)*sincY)^2;

   else
      y = 0;
   end
else
   Z = 2*pi*aca*sin(theta);
   if Z == 0
      besZ = 0.5;
   else
      besZ = besselj(1,Z)/Z;
   end
   
   if (adis == 1 & agr == 1)
      y = (sin(phi)^2+(cos(theta)*cos(phi))^2)* ...
         (besZ)^2;
   elseif (adis == 1 & agr == 2)
      y = ((1+cos(theta))/2*besZ)^2;
   elseif (adis == 2 & agr == 1)
      y = (sin(phi)*besZ)^2+ ...
         (cos(theta)*cos(phi)*(besselj(0,Z)-besZ)/(1-(Z/1.841)^2))^2;
   elseif (adis == 2 & agr == 2)
      y = ((1+cos(theta))/2)^2*((sin(phi)*besZ)^2+(cos(phi)*(besselj(0,Z)-besZ)/(1-(Z/1.841)^2))^2);
   else
      y = 0;
   end
end


function y = bsect(eh,left,right)
epsilon = 1e-4;
if eh == 3
   vl = fnbw(left);
   vr = fnbw(right);
else
   vl = hpbw(eh,left);
   vr = hpbw(eh,right);
end

if (vl*vr) > 0
   fprintf('exception in the program, quiting ... ');
end

while abs(left-right)>epsilon
   mid = (left+right)/2;
   if eh == 3
      vm = fnbw(mid);
   else
      vm = hpbw(eh,mid);
   end
      
   if (vm*vr)<0
      left = mid;
      vl = vm;
   else
      right = mid;
      vr = vm;
   end
end
y = (left+right)/2;


function y = hpbw(eh,theta)
global atype adis agr ara arb aca;

if atype == 1
   X = pi*ara*sin(theta);
   Y = pi*arb*sin(theta);
   if (adis == 1 & agr == 1)
      if (eh == 1)
         y = sin(Y)/Y;

      else
         y = cos(theta)*sin(X)/X;
      end
      
      y = abs(y) - sqrt(0.5);
      
   elseif (adis == 1 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*sin(Y)/Y;
      else
         y = (1+cos(theta))/2*sin(X)/X;
      end
      
      y = abs(y) - sqrt(0.5);
      
   elseif (adis == 2 & agr == 1)
      if (eh == 1)
         y = sin(Y)/Y;
         y = abs(y) - sqrt(0.5);
      else
         y = cos(theta)*cos(X)/(X^2-(pi/2)^2);
         y = abs(y) - (4/pi^2)/sqrt(2);
      end
   elseif (adis == 2 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*sin(Y)/Y;
         y = abs(y) - sqrt(0.5);
      else
         y = (1+cos(theta))/2*cos(X)/(X^2-(pi/2)^2);
         y = abs(y) - (4/pi^2)/sqrt(2);
      end
   end
   
elseif atype == 2
   Z = 2*pi*aca*sin(theta);
   if (adis == 1 & agr == 1)
      if (eh == 1)
         y = besselj(1,Z)/Z;
      else
         y = cos(theta)*besselj(1,Z)/Z;
      end
   elseif (adis == 1 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*besselj(1,Z)/Z;
      else
         y = (1+cos(theta))/2*besselj(1,Z)/Z;
      end
      
      
   elseif (adis == 2 & agr == 1)
      if (eh == 1)
         y = besselj(1,Z)/Z;
	      else
         y = cos(theta)*(besselj(0,Z)-besselj(1,Z)/Z)/(1-(Z/1.841)^2);
         
      end
   elseif (adis == 2 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*besselj(1,Z)/Z;
      else
         
         y = (1+cos(theta))/2*(besselj(0,Z)-besselj(1,Z)/Z)/(1-(Z/1.841)^2);
         
      end
      
   end
   
   y = abs(y) - 0.5/sqrt(2);   
   

end

function y = fnbw(theta)   
global atype adis agr ara arb aca;

if atype == 1
   X = pi*ara*sin(theta);
   y = cos(theta)*cos(X)/(X^2-(pi/2)^2);
elseif atype == 2
   Z = 2*pi*aca*sin(theta);
   y = cos(theta)*(besselj(0,Z)-besselj(1,Z)/Z)/(1-(Z/1.841)^2);
end


function y = amp(eh,theta)
global atype adis agr ara arb aca;

if atype == 1
   X = pi*ara*sin(theta);
   if X == 0
      sincX = 1;
   else
      sincX = sin(X)/X;
   end
   
   Y = pi*arb*sin(theta);
   if Y == 0
      sincY = 1;
   else
      sincY = sin(Y)/Y;
   end
   
   if (adis == 1 & agr == 1)
      if (eh == 1)
         y = sincY;

      else
         y = cos(theta)*sincX;
      end
      
%      y = abs(y) - sqrt(0.5);
      
   elseif (adis == 1 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*sincY;
      else
         y = (1+cos(theta))/2*sincX;
      end
      
      %y = abs(y) - sqrt(0.5);
      
   elseif (adis == 2 & agr == 1)
      if (eh == 1)
         y = sincY;
%         y = abs(y) - sqrt(0.5);
      else
         y = cos(theta)*cos(X)/(X^2-(pi/2)^2);
 %        y = abs(y) - (4/pi^2)/sqrt(2);
		end
      
   elseif (adis == 2 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*sincY;
%         y = abs(y) - sqrt(0.5);
      else
         y = (1+cos(theta))/2*cos(X)/(X^2-(pi/2)^2);
 %        y = abs(y) - (4/pi^2)/sqrt(2);
      end

   end

elseif atype == 2
   Z = 2*pi*aca*sin(theta);
   if Z == 0
      besZ = 0.5;
   else
      besZ = besselj(1,Z)/Z;
   end
   
   if (adis == 1 & agr == 1)
      if (eh == 1)
         y = besZ;
      else
         y = cos(theta)*besZ;
      end
   elseif (adis == 1 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*besZ;
      else
         y = (1+cos(theta))/2*besZ;
      end
      
      
   elseif (adis == 2 & agr == 1)
      if (eh == 1)
         y = besZ;
         
      else
         y = cos(theta)*(besselj(0,Z)-besZ)/(1-(Z/1.841)^2);
      end
      
   elseif (adis == 2 & agr == 2)
      if (eh == 1)
         y = (1+cos(theta))/2*besZ;
         
      else
         y = (1+cos(theta))/2*(besselj(0,Z)-besZ)/(1-(Z/1.841)^2);
      end

   end
%   y = abs(y) - 0.5/sqrt(2);   
end
y = abs(y);



function y = func3d(z)
y = z;
for i = 1:size(z,1)
   for j = 1:size(z,2)
      if real(z(i,j)) == 0
         phi = 0;
      else
         phi = atan(imag(z(i,j))/real(z(i,j)));
      end
      
      theta = abs(z(i,j))*pi/2;
      y(i,j) = abs(U(theta,phi));
   end
end





function hpol = polar_dB(theta,rho,rmin,rmax,rticks,line_style) 

% Font size, font style and line width parameters
font_size  = 16;
font_name  = 'Times';
line_width = 1.5;

if nargin < 5
	error('Requires 5 or 6 input arguments.')
elseif nargin == 5 
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 1
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   font_size, ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
  % v returns the axis limits
  % changed the following line to let the y limits become negative
	hhh=plot([0 max(theta(:))],[min(rho(:)) max(rho(:))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);

% check radial limits (rticks)

 	if rticks > 5   % see if we can reduce the number
 		if rem(rticks,2) == 0
 			rticks = rticks/2;
 		elseif rem(rticks,3) == 0
 			rticks = rticks/3;
 		end
 	end

% define a circle
	th = 0:pi/50:2*pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:(length(th)-1)/4:length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);

	rinc = (rmax-rmin)/rticks;

% label r
  % change the following line so that the unit circle is not multiplied
  % by a negative number.  Ditto for the text locations.
	for i=(rmin+rinc):rinc:rmax
                is = i - rmin;
		plot(xunit*is,yunit*is,'-','color',tc,'linewidth',0.5);
		text(0,is+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
	end
% plot spokes
	th = (1:6)*2*pi/12;
	cst = cos(th); snt = sin(th);
	cs = [-cst; cst];
	sn = [-snt; snt];
	plot((rmax-rmin)*cs,(rmax-rmin)*sn,'-','color',tc,'linewidth',0.5);

% plot the ticks
	george=(rmax-rmin)/30; % Length of the ticks
        th2 = (0:36)*2*pi/72;
        cst2 = cos(th2); snt2 = sin(th2);
	cs2 = [(rmax-rmin-george)*cst2; (rmax-rmin)*cst2];
	sn2 = [(rmax-rmin-george)*snt2; (rmax-rmin)*snt2];
	plot(cs2,sn2,'-','color',tc,'linewidth',0.15); % 0.5
        plot(-cs2,-sn2,'-','color',tc,'linewidth',0.15); % 0.5


% annotate spokes in degrees
  % Changed the next line to make the spokes long enough
	rt = 1.1*(rmax-rmin);
	for i = 1:max(size(th))
		text(rt*cst(i),rt*snt(i),int2str(abs(i*30-90)),'horizontalalignment','center' );
		if i == max(size(th))
			loc = int2str(90);
		elseif i*30+90<=180
			loc = int2str(i*30+90);
                else
                        loc = int2str(180-(i*30+90-180));  
		end
		text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center' );
	end
% set viewto 2-D
	view(0,90);

% set axis limits
  % Changed the next line to scale things properly
	axis((rmax-rmin)*[-1 1 -1.1 1.1]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
  % changed the next line so negative rho are not plotted on the other side
  
  for i = 1:length(rho)
    if (rho(i) > rmin)
      if theta(i)*180/pi >=0 & theta(i)*180/pi <=90
          xx(i) = (rho(i)-rmin)*cos(pi/2-theta(i));
          yy(i) = (rho(i)-rmin)*sin(pi/2-theta(i));
      elseif theta(i)*180/pi >=90
          xx(i) = (rho(i)-rmin)*cos(-theta(i)+pi/2);
          yy(i) = (rho(i)-rmin)*sin(-theta(i)+pi/2);
      elseif theta(i)*180/pi < 0 
          xx(i) = (rho(i)-rmin)*cos(abs(theta(i))+pi/2);
          yy(i) = (rho(i)-rmin)*sin(abs(theta(i))+pi/2);
      end
    else
      xx(i) = 0;
      yy(i) = 0;
    end
  end

% plot data on top of grid
if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style);
end
if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end
