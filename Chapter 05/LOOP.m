%***********************************************************************
% LOOP.m
%**********************************************************************
% This is a MATLAB based program computes the:
%	I.   Maximum directivity (dimensionless and in dB)
%	II.	 Radiation resistance
%   III. Directivity pattern (in dB) in polar form
%	IV.  Normalized elevation far-field pattern (E-phi, in dB)
%
%   and also plots the:
% 	I.   Directivity polar pattern (in dB)
%	II.  Normalized far-field polar amplitude pattern (E-phi, in dB) 
%
% for a circular loop (with constant current).  The loop is radiating
% into free space.
%
% The directivity and radiation resistance are calculated using
% the trailing edge method in increments of 1 degree in theta.
%
%	**Input parameters:
%	1.	A: Loop radius (in wavelengths)
%
%	**Note:
%	The far-zone electric field component E-phi exists for
%	0 < theta < 180 and 0 < phi < 360.
%-----------------------------------------------------------------
% Converted from Fortran to Matlab 3/2002 by Kelly O'Dell.
% Modified by Marios Gkatzianas
%-----------------------------------------------------------------
%
clear all;
close all;
format long;

warning off;

%---Choice of output---

fprintf('Output device option \n\tOption (1): Screen\n\tOption (2): File \n');
ERR = 1;
while(ERR ~= 0)
   DEVICE = input('\nOutput device = ','s');
   DEVICE = str2num(DEVICE);
   if(DEVICE == 1)
      ERR = 0;
   elseif(DEVICE == 2)
      FILNAM = input('Input the desired output filename: ','s');
      ERR = 0;
   else
      error('Outputting device number should be either 1 or 2\n');
   end
end


%---Definition of constants and initialization---

PI = 4.0*atan(1.0);
E = 120.0*PI;
THETA = PI/180.0;
UMAX = 0.0;
PRAD = 0.0;
TOL = 1.0E-5;

%---Input the radius of the loop---

A = input('\nRadius of loop in wavelengths = ','s');
A = str2num(A);
%***Insert input data error loop*****

%---Main program------------------------------------

I = 1;
while(I <= 180)
   XI = I*PI/180.0;
   X = 2.0*PI*A*sin(XI);
   if(abs(X) < TOL)
      F = 0.0;
   else
      F = besselj(1,X);
   end
   U = A^2*(2.0*PI)^2/8.0*E*F^2;
   if(U > UMAX)
      UMAX = U;
   end
   UA = U*sin(XI)*THETA*2*PI;
   PRAD = PRAD+UA;
   I = I+1;
end
D = (4.0*PI*UMAX)/PRAD;
DDB = 10.0*log10(D);
RR = 2.0*PRAD;

%---Calculation of elevation far-field patterns in 1 degree increments---
fid = fopen('ElevPat.dat','w');
fprintf(fid,'\tLoop\n\n\tTheta\t\tH (dB)\n');
fprintf(fid,'\t-----\t\t------\n');
T = zeros(180,1);
HT = zeros(180,1);
HdB = zeros(180,1);

x = 1;
while(x<=180)
   T(x) = x-0.99;
   Y = 2*PI*A*sin(T(x)*THETA);
   HT(x) = besselj(1,Y);
   x = x+1;
end
HT = abs(HT);
HTmax = max(HT);
HdB = 20*log10(abs(HT)/HTmax);

x = 1;
while(x<=180)
   fprintf(fid,'\n %5.0f %12.4f',T(x),HdB(x));
   x = x+1;
end

fclose(fid);


%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end

x=linspace(0,4*pi*A,500);
dx=x(2)-x(1);
den=sum(besselj(2,x)*dx);
theta=linspace(0,pi,300);

Dth=4*pi*A/den*(besselj(1,2*pi*A*sin(theta))).^2;
Dth(1)=0;
Dth_db=10*log10(Dth);


%---Echo input parameters and output computed parameters---
fprintf(fid,'\nLOOP:\n-----');
fprintf(fid,'\n\nInput parameters:\n-----------------');
fprintf(fid,'\nRadius of loop in wavelengths = %6.4f',A);
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nDirectivity (dimensionless) = %6.4f',D);
fprintf(fid,'\nDirectivity (dB) \t= %6.4f',DDB);
fprintf(fid,'\nRadiation resistance (Ohms) = %9.4f',RR);
fprintf(fid,'\n\n***NOTE:\nThe normalized elevation pattern is stored\n');
fprintf(fid,'in an output file called ..........ElevPat.dat\n\n');

if(DEVICE == 2)
   fclose(fid);
end

%---Plot elevation far field pattern------
% plot(T,-abs(HdB),'b');
% axis([0 180 -60 0]);
% grid on;
% xlabel('Theta (degrees)');
% ylabel('Amplitude (dB)');
% legend(['a = ',num2str(A),' \lambda'],0);
% title('Loop Far-Field Elevation Pattern');

T=T'; HdB=HdB';

T=[T T+180];
Hdb=[HdB fliplr(HdB)];

% Figure 1
% ********
polar_dB(T,-abs(Hdb),-60,max(Hdb),4);
title('Elevation plane normalized amplitude pattern (dB)','fontsize',16);

% Figure 2
% ********
figure(2);
polar_dB([theta theta+pi]*180/pi,[Dth_db fliplr(Dth_db)],-60,max(Dth_db),4);
title('Elevation plane directivity pattern (dB)','fontsize',16);

%---End of program-----------------------------------------
