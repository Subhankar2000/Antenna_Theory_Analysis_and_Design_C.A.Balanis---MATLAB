% ************************************************************************************
% Polar
%**************************************************************************
% This interactive MATLAB program is used to plot two-dimensional (2-D):
%
% A.  Polar (-180 degrees < PATTERN < +180 degrees)
% or
% B.  Semipolar (-90 degrees < PATTERN < +90 degrees)
%
% To accomplish this, PATTERN calls Matlab functions:
%
% C.  Polar_dB.m or Semipolar_dB.m
% D.  f.m
%
% The user MUST specify the function, which represents the LINEAR FIELD pattern 
% of the function whose pattern is to be plotted, in the Matlab function
%
% f.m
%
% The program can plot any of the following:
%
% 1.  Field pattern (linear scale)
% 2.  Power pattern (linear scale)
% 3.  Power pattern (in dB)
%
% Program was written by:
% Bo Yang, AHE Lab 
% Department of Electrical Engineering
% Arizona State University
% February 2004


function Patt2D

clear all;
close all;

% Data Entry
fprintf('\n-------------------------------------------------------');
fprintf('\nThis program calls function Polar_dB.m or Semipolar_dB.m');
fprintf('\nto plot 2D Pattern in Polar Coordinates');
fprintf('\n-------------------------------------------------------\n');
fprintf('\n   *** NOTICE: Please edit m-file "f.m" to input field pattern expressions!\n\n');
fprintf('   Make sure to use dot operations(".*","./" and ".^") in the expression,\n');
fprintf('   instead of regular operations("*","/" and "^")!\n\n\n');;

% Plot Range
ERR = 1;
while(ERR ~= 0)
   fprintf('\nPlease specify the range of theta in the form "[theta1,theta2]"\n');
   fprintf('(-180 <= theta1 <= theta2 <= 180 (degrees), for example: [-180,180])\n');
   fprintf('-----------------------------------------------------------------\n');
   thrg = str2num(input('Input: ','s'));
   s=size(thrg);
   if (s(1) ~= 1) | (s(2) ~= 2)
      fprintf('\n   *** ERROR: Please specify the range in the form "[theta1,theta2]")!\n\n','r');
   elseif (thrg(1) < -180) | (thrg(1) >= thrg(2)) | (thrg(2) > 180)
      fprintf('\n   *** ERROR: Theta must satisfys: -180 <= theta1 <= theta2 <= 180 (degrees)")!\n\n');
   else
      ERR = 0;
   end
end

% Angle Increment
ERR = 1;
while(ERR ~= 0)
   fprintf('\nPlease specify the angle increment (in degrees)"\n');
   fprintf('--------------------------------------------------\n');
   del_th = str2num(input('Input: ','s'));
   s=size(del_th);
   if (s(1) ~= 1) | (s(2) ~= 1)
      fprintf('\n   *** ERROR: Please specify the angle increment")!\n\n');
   elseif (del_th < 0) | (del_th >= 180)
      fprintf('\n   *** ERROR: del_th must satisfy: 0 <= del_th <= 180 (degrees)")!\n\n');
   else
      ERR = 0;
   end
end

thrg(2)=round((thrg(2)-thrg(1))/del_th)*del_th+thrg(1);
theta=((thrg(1):del_th:thrg(2)));

% Select plotting scale
ERR = 1;
rmin = 0;
while(ERR ~= 0)
    fprintf('\nPlease select the display scale.\n');
    fprintf('----------------------------------\n');
    fprintf('(1) Field Pattern ( linear scale ).\n');
    fprintf('(2) Power Pattern ( linear scale ).\n');    
    fprintf('(3) Power Pattern ( in dB ).\n\n');
    scl = str2num(input('Input: ','s'));
    if (scl == 3)
        ERR1=1;
        while(ERR1 ~= 0)
           fprintf('\nPlease enter the minimum limit (in dB)of the plot (e.g., -40 dB)\n');
           fprintf('------------------------------------------------------------------\n');           
           rmin = str2num(input('Input:','s'));
           if rmin < 0
               ERR1 = 0;
           else
               fprintf('\n   *** ERROR: Please enter a negtive number!\n\n');
           end
        end
        ERR = 0;
    elseif (scl == 1) | (scl == 2)
        ERR = 0;
    end
end

% Select plotting Scheme
ERR = 1;
while(ERR ~= 0)
    fprintf('\nPlease select the plotting scheme.\n');
    fprintf('------------------------------------\n');
    fprintf('(1) Polar.\n');
    fprintf('(2) Semipolar.\n');    
    pltsch = str2num(input('Input: ','s'));
    if (pltsch == 1) | (pltsch == 2)
        ERR = 0;
    end
end

% Calculating
fprintf('Calulating ... \n\n\n');
r=abs(f(theta*pi/180));

% Scaling
r=r/max(r);
if scl < 3
   r=r.^scl;
else
   r=20*log10(r);
end

% Plot
if pltsch == 1
    polar_dB(theta,r,rmin,max(r),4,'-');
else
    idx=find((theta >= -90)&(theta <= 90)); 
    semipolar_dB(theta(idx),r(idx),rmin,max(r),4,'-');
end
