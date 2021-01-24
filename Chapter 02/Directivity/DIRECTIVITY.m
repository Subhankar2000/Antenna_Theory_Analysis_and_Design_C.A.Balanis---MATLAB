%*******************************************************************************************
% DIRECTIVITY.m
%*************************************************************************
% This is a MATLAB based program that calculates the: 
%	I. 	radiated power 
%	II.	maximum directivity (dimensionless and in dB) 
% of any antenna.  The maximum directivity is calculated using the trailing 
% edge method in increments of 1 degree in theta and 1 degree in phi.
% 
% Input parameters:
%	1.  TL, TU: Lower and upper limits in theta (in degrees)
%	2.  PL, PU: Lower and upper limits in phi (in degrees)
%	3.  U(theta,phi): The radiation intensity function (in file 'U.m')
%
%	**NOTE: The radiation intensity function F must be provided for a given
%				antenna and should be inserted into the subroutine U.
%
%	**EXAMPLE: If the antenna is radiating only in the upper hemisphere,
%		the input parameters are TL = 0, TU = 90, PL = 0, PU = 360.
%---------------------------------------------------------------------------------------
% Converted from fortran to  MATLAB 3/2002 by Kelly O'Dell
% Modified by Marios Gkatzianas
% ---------------------------------------------------------------------------------------

close all;
clear all;
format long;

fprintf('!!!WARNING: Make sure you define the radiation intensity in the file U.m !!! \n\n');

s=which('U.m');
if isempty(s),
   disp(['File U.m not found in current directory. Make sure that the radiation intensity is' ...
          ' defined in that file!!']);         
end;   

device=3;

while ((device~=1)&(device~=2)),
   fprintf('Output device option \n\tOption (1): Screen\n\tOption (2): File \n\n');
   device=input('Output device=');
   
   if device==2,
      filename=input('Input the desired output filename: ','s');
   elseif device~=1,
      fprintf('\nOutputting device number should be either 1 or 2\n');
   end;   
end;

TL=[]; TU=[]; PL=[]; PU=[];

while isempty(TL),
   TL = input('\nThe lower bound of theta in degrees = ');
end;   

while isempty(TU),
   TU = input('The upper bound of theta in degrees = ');
end;

while isempty(PL),
   PL = input('The lower bound of phi in degrees = ');
end;

while isempty(PU),
   PU = input('The upper bound of phi in degrees = ');
end;

theta=(TL:TU)*pi/180;  % convert into radians
phi=(PL:PU)*pi/180;

dth=theta(2)-theta(1);
dph=phi(2)-phi(1);

[THETA,PHI]=meshgrid(theta,phi);

x=U(THETA,PHI);

Prad=sum(sum(x.*sin(THETA)*dth*dph));
D=4*pi*max(max(x))/Prad;
D_db=10*log10(D);

if (device==2),
   fid=fopen(filename,'wt');
else
   fid=2;   % standard output   
end;

fprintf(fid,'\nInput parameters:\n-----------------');
fprintf(fid,'\nThe lower bound of theta in degrees = %3.0f',TL);   
fprintf(fid,'\nThe upper bound of theta in degrees = %3.0f',TU);
fprintf(fid,'\nThe lower bound of phi in degrees = %3.0f',PL);
fprintf(fid,'\nThe upper bound of phi in degrees = %3.0f',PU);
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nRadiated power (watts) = %6.4f',Prad);
fprintf(fid,'\nDirectivity (dimensionless) = %6.4f',D);
fprintf(fid,'\nDirectivity (dB) \t= %6.4f',D_db);
fprintf(fid,'\n\n');

if(device == 2)
   fclose(fid);
end

