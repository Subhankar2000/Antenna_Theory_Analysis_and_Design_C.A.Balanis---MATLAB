%***********************************************************************
%     Design.m
%***********************************************************************
%     THIS PROGRAM DESIGNS AN OPTIMUM GAIN PYRAMIDAL HORN BASED ON THE
%     FORMULATION OF SECTION 13.4.3, EQUATIONS (13-55) - (13-58b).
%
%     THE PROGRAM COMPUTES THE:
%
%     I.   HORN DIMENSIONS A1, B1, RHOE, RHOH, PE, PH (IN CM)
%     II.  HORN FLARE ANGLES PSIE (E-PLANE) AND PSIH (H-PLANE) 
%          (IN DEGREES)
%
%       ** INPUT PARAMETERS:
%
%       1.  DESIRED GAIN Go (IN dB)
%       2.  FREQUENCY (IN GHz)
%       3.  FEED RECTANGULAR WAVEGUIDE DIMENSIONS: A, B (IN CM)
%
%       ** NOTE: REFER TO FIGURE 13.18 FOR THE GEOMETRY
%
%******************************************************************
%       Written by: Steven Parker, Arizona State University
%
%******************************************************************       
%-------------------------------------------------------------
%   Ask User for Input Parameters
%-------------------------------------------------------------

clear;
global G a b Lambda f;
fprintf('\nOUTPUT DEVICE OPTION\n');
fprintf('\tOPTION (1): SCREEN\n');
fprintf('\tOPTION (2): OUTPUT FILE\n');
output=input('OUTPUT DEVICE:  ');
fprintf('\n---------------------------------------');
fprintf('\nSpecify the following input parameters:\n');
fprintf('---------------------------------------');

Go=input('\nDESIRED GAIN OF THE HORN IN dB: Go(dB)= ');
fo=input('FREQUENCY OF OPERATION IN GHz: fo(GHz)= ');
a=input('HORN DIMENSION A IN CM: a(cm)= ');
b=input('HORN DIMENSION B IN CM: b(cm)= ');


%----------------------------------------------
%   Define Variables: G, f, Lambda, a & b
%----------------------------------------------

f=fo*1e9;
G=10^(Go/10);
Lambda=((3e8)/f)*100;

%----------------------------------------------
%   Compute Chi by solving Eq.(13-56)
%   Begin with initial Chi value given by
%   Eq.(13-57)
%----------------------------------------------

X1=G/(2*pi*sqrt(2*pi));         %Eq.(13-57)
Chi=fzero(@Subroutine,X1);     %Use subroutine

%----------------------------------------------
%   Compute RHOe & RHOh: Units (cm) 
%   Eq.(13-56a) & (13-56b)
%----------------------------------------------

RHOe=Chi*Lambda;
RHOh=((G^2)/(8*pi^3))*(1/Chi)*Lambda;

%----------------------------------------------
%   Compute a1 & b1: Units (cm)
%   Eq.(13-56a) & (13-56b)
%----------------------------------------------

a1=(G/(2*pi))*sqrt(3/(2*pi*Chi))*Lambda;
b1=sqrt(2*Chi)*Lambda;

%----------------------------------------------
%   Compute Pe & Ph: Units (cm)
%   Eq.(13-49a) & (13-49b)
%----------------------------------------------

Pe=(b1-b)*sqrt((RHOe/b1)^2-1/4);
Ph=(a1-a)*sqrt((RHOh/a1)^2-1/4);

%----------------------------------------------
%   Compute PSIe & PSIh: Units (degrees)
%   Refer to Figure 13.18 for Geometry
%----------------------------------------------

PSIe=asin(b1/(2*RHOe))*(180/pi);
PSIh=asin(a1/(2*RHOh))*(180/pi);

%---------------------------------------------------
%   Display Optimum Gain Pyramidal Horn Parameters
%---------------------------------------------------

switch output
case 1
fprintf('\n---------------------------------------------');
fprintf('\nDESIGNED PARAMETERS FOR THE OPTIMUM GAIN HORN');
fprintf('\n---------------------------------------------');
fprintf('\na1\t\t=\t%f',a1);
fprintf('\tcm');
fprintf('\nb1\t\t=\t%f',b1);
fprintf('\tcm');
fprintf('\nRHOe\t=\t%f',RHOe);
fprintf('\tcm');
fprintf('\nRHOh\t=\t%f',RHOh);
fprintf('\tcm');
fprintf('\nPe\t\t=\t%f',Pe);
fprintf('\tcm');
fprintf('\nPh\t\t=\t%f',Ph);
fprintf('\tcm');
fprintf('\nPSIe\t=\t%f',PSIe);
fprintf('\tDeg');
fprintf('\nPSIh\t=\t%f',PSIh);
fprintf('\tDeg');

case 2

fid = fopen(input('Specify file name: ','s'),'wt');
fprintf(fid,'*********PROGRAM OUTPUT*********\n');
fprintf(fid,'INPUT PARAMETERS\n');
fprintf(fid,'----------------\n');
fprintf(fid,'DESIRED GAIN OF THE HORN\t =\t%f',Go);
fprintf(fid,'\tdB');
fprintf(fid,'\nFREQUENCY OF OPERATION\t\t =\t%f',fo);
fprintf(fid,'\tGHz');
fprintf(fid,'\nFEED WAVEGUIDE DIMENSION A\t =\t%f', a);
fprintf(fid,'\tcm');
fprintf(fid,'\nFEED WAVEGUIDE DIMENSION B\t =\t%f', b);
fprintf(fid,'\tcm');
fprintf(fid,'\n\nDESIGNED PARAMETERS FOR THE OPTIMUM GAIN HORN');
fprintf(fid,'\n---------------------------------------------\n');

fprintf(fid,'a1\t=\t%f',a1);
fprintf(fid,'\tcm');
fprintf(fid,'\nb1\t=\t%f',b1);
fprintf(fid,'\tcm');
fprintf(fid,'\nRHOe\t=\t%f',RHOe);
fprintf(fid,'\tcm');
fprintf(fid,'\nRHOh\t=\t%f',RHOh);
fprintf(fid,'\tcm');
fprintf(fid,'\nPe\t=\t%f',Pe);
fprintf(fid,'\tcm');
fprintf(fid,'\nPh\t=\t%f',Ph);
fprintf(fid,'\tcm');
fprintf(fid,'\nPSIe\t=\t%f',PSIe);
fprintf(fid,'\tDeg');
fprintf(fid,'\nPSIh\t=\t%f',PSIh);
fprintf(fid,'\tDeg');



status = fclose(fid);

end