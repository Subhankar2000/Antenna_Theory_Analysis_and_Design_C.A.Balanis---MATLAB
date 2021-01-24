
% ****************************************************************************
% DIPOLE.m
%*************************************************************************
% This is a MATLAB based program that computes the: 
%	I.	 Maximum directivity (dimensionless and in dB)
%	II.	 Radiation resistance (Rr)
%   III. Input resistance (Rin)
%	IV.  Reactance relative to current maximum (Xm)
%   V.   Input reactance (Xin)
%	VI.  Normalized current distribution
%   VII. Directivity pattern (in dB) in polar form
%   VIII.Normalized far-field amplitude pattern (E-theta, in dB) in polar form
% for a symmetrical dipole of finite length.  The dipole is radiating
% in free space.
%
% The directivity, resistances and resistances are calculated using the trailing 
% edge method in increments of 1 degree in theta.
%
% **Input parameters
% 1.	L:	Dipole length (in wavelengths)
% 2.    a:  Dipole radius (in wavelengths)
%
% **Note:
% The far zone electrif field component, E-theta, exists for
% 0 < theta < 180 and 0 < phi < 360.
%------------------------------------------------------------------------
% Converted from Fortran to Matlab by Kelly O'Dell 3/2002
% Modified by Marios Gkatzianas
%------------------------------------------------------------------------
%
function []=dipole;

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
TOL = 1.0E-6;

%---Input the length of the dipole---

L = input('\nLength of dipole in wavelengths = ','s');
L = str2num(L);
%***Insert input data error loop***

r=input('Radius of dipole in wavelengths = ');


%---Main program----------------------

A = L*PI;
I = 1;
while(I <= 180)
   XI = I*PI/180.0;
   if(XI ~= PI)
      U = ((cos(A*cos(XI))-cos(A))/sin(XI))^2*(E/(8.0*PI^2));
      if(U > UMAX)
         UMAX = U;
      end
   end
   UA = U*sin(XI)*THETA*2.0*PI;
   PRAD = PRAD+UA;
   I = I+1;
end
D = (4.0*PI*UMAX)/PRAD;
DDB = 10.0*log10(D);
RR = 2.0*PRAD;
if(A ~= PI)
   RIN = RR/(sin(A))^2;
end


%---Calculation of elevation far-field patterns in 1 degree increments---
fid = fopen('ElevPat.dat','w');
fprintf(fid,'\tDipole\n\n\tTheta\t\tE (dB)\n');
fprintf(fid,'\t----\t\t------');
T = zeros(180,1);
ET = zeros(180,1);
EdB = zeros(180,1);

x = 1;
while(x<=180)
   T(x) = x-0.99;
   ET(x) = (cos(PI*L*cos(T(x)*THETA))-cos(PI*L))/sin(T(x)*THETA);
   x = x+1;
end
ET = abs(ET);
ETmax = max(abs(ET));
EdB = 20*log10(abs(ET)/ETmax);

x = 1;
while(x<=180)
   fprintf(fid,'\n %5.4f %12.4f',T(x),EdB(x));
   x = x+1;
end

fclose(fid);


n=120*pi;
k=2*pi;

if exist('cosint')~=2,
   disp(' ');
   disp('Symbolic toolbox is not installed. Switching to numerical computation of sine and cosine integrals.');
   Xm=30*(2*si(k*L)+cos(k*L)*(2*si(k*L)-si(2*k*L))- ...
          sin(k*L)*(2*ci(k*L)-ci(2*k*L)-ci(2*k*r^2/L)));
   Xin=Xm/(sin(k*L/2))^2;
      
elseif exist('cosint')==2,
   Xm=30*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))- ...
          sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint(2*k*r^2/L)));
   Xin=Xm/(sin(k*L/2))^2;
end;   



%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end

%---Echo input parameters and output computed parameters---
fprintf(fid,'\nDIPOLE:\n-------');
fprintf(fid,'\n\nInput parameters:\n-----------------');
fprintf(fid,'\nLength of dipole in wavelengths = %6.4f',L);
fprintf(fid,'\nRadius of dipole in wavelengths = %6.7f',r);
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nDirectivity (dimensionless) = %6.4f',D);
fprintf(fid,'\nDirectivity (dB) \t= %6.4f\n',DDB);
fprintf(fid,'\nRadiation resistance based on current maximum (Ohms) = %10.4f',RR);
fprintf(fid,'\nReactance based on current maximum (Ohms) = %10.4f\n',Xm);

if(abs(sin(A)) < TOL)
   fprintf(fid,'\nInput resistance = INFINITY');
   fprintf(fid,'\nInput reactance = INFINITY\n\n');   
else
   fprintf(fid,'\nInput resistance (Ohms) = %10.4f',RIN);
 % fprintf(fid,'\nInput reactance based on current maximum (Ohms) = %10.4f',Xm);
   fprintf(fid,'\nInput reactance (Ohms) = %10.4f',Xin);
   fprintf(fid,'\n\n***NOTE:\nThe normalized elevation pattern is stored\n');
   fprintf(fid,'in an output file called ..........ElevPat.dat\n\n');
end


if(DEVICE == 2)
   fclose(fid);
end

%---Plot elevation far field pattern------
% plot(T,EdB,'b');
% axis([0 180 -60 0]);
% grid on;
% xlabel('Theta (degrees)');
% ylabel('Amplitude (dB)');
% legend(['L = ',num2str(L),' \lambda'],0);
% title('Dipole Far-Field Elevation Pattern');


% Figure 1
% ********
z=linspace(-L/2,L/2,500);
k=2*pi;
I=sin(k*(L/2-abs(z)));
plot(z,abs(I));
xlabel('z^{\prime}/\lambda','fontsize',12);
ylabel('Normalized current distribution','fontsize',12);


% Figure 2
% ********
figure(2);
T=T'; EdB=EdB';
EdB=[EdB fliplr(EdB)];
T=[T T+180];

polar_dB(T,EdB,-60,0,4);
title('Elevation plane normalized amplitude pattern (dB)','fontsize',16);

% Figure 3
% ********
figure(3);
theta=linspace(0,2*pi,300);

Eth=(cos(k*L/2*cos(theta))-cos(k*L/2))./sin(theta);
Dth=4*pi*120*pi/(8*pi^2)*Eth.^2/PRAD;
Dth_db=10*log10(Dth);
Dth_db(Dth_db<=-60)=-60;

polar_dB(theta*180/pi,Dth_db,-60,max(Dth_db),4);
title('Elevation plane directivity pattern (dB)','fontsize',16);

%---End program----------------------------------------------


function [y]=si(x);

v=linspace(0,x/pi,500);
dv=v(2)-v(1);

y=pi*sum(sinc(v)*dv);


function [y]=ci(x);

v=linspace(0,x/(2*pi),500);
dv=v(2)-v(1);

y1=2*pi*sum(sinc(v).*sin(pi*v)*dv);

y=.5772+log(x)-y1;



