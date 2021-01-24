% ********************************************************************
%    Beamforming
%*********************************************************************
%
%    THIS IS A MATLAB BASED PROGRAM THAT COMPUTES THE WEIGHTS AND
%    BEAMFORMED PATTERN OF:
%
%     I.   LINEAR ANTENNA
%     II.  RECTANGULAR PLANAR ANTENNA
%
%     THE LINEAR ARRAY HAS M ELEMENTS PLACED EQUIDISTANTLY ALONG
%     THE X-AXIS. 
%
%     THE RECTANGULAR PLANAR ARRAY HAS M x N ELEMENTS PLACED EQUIDISTANTLY
%     ALONG THE X AND Y AXES.
%     
%     OPTION I.  LINEAR ARRAY
%
%       ** INPUT PARAMETERS:
%
%          1. NUMBER OF ELEMENTS
%          2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%          3. AMPLITUDE AND DIRECTION (THETA)OF SIGNAL-OF-INTEREST(SOI)
%          4. AMPLITUDE AND DIRECTION (THETA)OF SIGNALS-NOT-OF-INTEREST
%             (SNOIs)
%          5. NOISE MEAN AND VARIANCE
%          6. MU OF LMS ALGORITHM
%
%       ** PROGRAM OUTPUT:
%
%          1. ARRAY WEIGHTS (AMPLITUDE AND PHASE)
%          2. BEAMFORMED PATTERN 
%
%  OPTION II.  RECTANGULAR PLANAR ARRAY
%
%       ** INPUT PARAMETERS:
%
%          1. NUMBER OF ARRAY ELEMENTS IN X-DIRECTION
%          2. NUMBER OF ARRAY ELEMENTS IN Y-DIRECTION
%          3. SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS)
%          4. SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS)
%          5. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNAL-OF-INTEREST
%             (SOI)
%          6. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNALS-NOT-OF-
%             -INTEREST (SNOIs) 
%          7. NOISE MEAN AND VARIANCE
%          8. MU OF LMS ALGORITHM
%
%       ** PROGRAM OUTPUT:
%
%          1. ARRAY WEIGHTS (AMPLITUDE AND PHASE)
%          2. BEAMFORMED PATTERN 
%
%     NOTES: ALL THE INPUT PARAMETERS ARE IN TERMS OF THE WAVELENGTH.
%************************************************************************
%     Written by: Salvatore Bellofiore, 2002
%     Revised by: Zhiyong Huang, Arizona State University,  10/2/2004
%************************************************************************

% Clean-up variables, figures
clear all;
close all;

fprintf('\n *** NOTICE: PLEASE CONTACT DR C.A BALANIS IF YOU HAVE ANY QUESTION OR COMMENT. THANKS.\n\n');

%---Choice of array, linear or planar---------------------------------------------
fprintf('Smart antenna structure option: \n\tOption (1): Linear\n\tOption (2): Planar \n');
ERR = 1;
while(ERR ~= 0)
   array = str2num(input('\n Array structure = ','s'));
   if(array == 1)
       fprintf('\n Linear Array is chosen.\n');
       ERR = 0;
       beamforming_linear;
   elseif(array == 2)
       fprintf('\n Planar Array is chosen.\n');
       ERR = 0;
       beamforming_planar;
   else
      fprintf('\n Antenna structure number should be either 1 or 2\n');
   end
end