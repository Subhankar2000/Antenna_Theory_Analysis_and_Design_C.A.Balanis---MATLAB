%****************************************************************************************
% TEE.m
%*********************************************************************
% This is a MATLAB based program to analyze and design a T-Match. 
% 
% The T-match connection (which is a particular case of the folded dipole) is employed as 
% a broadband matching technique to feed a dipole from a mismatched, but most likely 
% unbalanced, transmission line (like a twin-lead line).  The inner, shorter length (and 
% usually smaller radius) dipole is connected to the transmission line at its center and 
% is tapped at two locations to the longer (outer) dipole, each tap point a distance l'/2 
% from its center.  The two dipoles are separated from each other by a small distance, s.  
% Two capacitors (of equal capacitance), whose purpose is to resonate the circuit and 
% balance the line, are placed at the feed points of the inner dipole.  (See Figure 9.20 
% from the book Antenna Theory by Balanis.)   
%
% This MATLAB program for a T-match connection is designed to calculate the values of
% l' and C, which are, respectively, the length of the inner dipole (fed by the
% transmission line) and the capacitance of each of the two capacitors at the junction 
% of the transmission line and the inner dipole.  
%
% The program assumes the geometry of the dipoles is fixed (other than the length of the
% inner dipole).  The user is prompted to select from one of two options: 
% EITHER a dipole antenna whose length is known 
% OR an arbitrary antenna whose input impedance is known.
% For either case the user is prompted to input the values of the two radii (whether the 
% two dipoles or the antenna with T-match connection) and the separation between the two 
% elements.  For the two-dipole case, the user inputs the length of the outer dipole; for 
% the antenna-T-match arrangement, the user inputs the impedance of the antenna (its real 
% and imaginary parts separately).  
%
% In addition, the user inputs the characteristic impedance (which is assumed to be real-
% valued) of the feed transmission line; this value is denoted Zc.  The user must also 
% input the frequency (in megahertz) at which the dipole is to resonate.
%
%
%-------------------------------WARNING!!-------------------------------------------------
% THIS PROGRAM WILL FAIL TO YIELD ANY SATISFACTORY ANSWER IF THE TRANSMISSION LINE 
% CHARACTERISTIC IMPEDANCE, Zc, IS TOO LARGE.  IN PARTICULAR, IF Zc EXCEEDS Zo/g2, THE 
% CALCULATED OUTPUT LENGTH L' WILL BE A COMPLEX NUMBER; THE VARIABLES Zo AND g2 ARE 
% CALCULATED WITHIN THE PROGRAM, AND THEIR VALUES DEPEND ON THE GEOMETRY SPECIFIED BY THE 
% USER (i.e., ON THE INPUT VALUES a, a', s, AND L).
%
% FURTHERMORE, IF Zc APPROACHES (FROM BELOW) TOO CLOSELY THE CRITICAL VALUE OF Zo/g2, THE 
% OUTPUT VALUE OF L' MAY BE REAL, BUT IT WILL BE NEGATIVE.
%
% THE CENTER-TO-CENTER SEPARATION, s, OF THE INNER AND OUTER DIPOLES MUST BE SUFFICIENTLY
% LARGE TO SATISFY THE CONSTRAINT s^2 > (a^2 + (a')^2).  THIS CONSTRAINT SHOULD NEVER POSE
% A PROBLEM FOR REALISTIC GEOMETRIES, BUT MUST BE HEEDED TO PREVENT THE USER FROM ENTERING
% TOO SMALL A VALUE OF s.
%-----------------------------------------------------------------------------------------
%
%
% The program computes and outputs the following values:
% 1. Effective radius of the dipole-dipole configuration, ae
% 2. Dipole input impedance, Za
% 3. Characteristic impedance of the dipole pair, Zo
% 4. Transmission line mode impedance, Zt
% 5. Total input resistance, Rin
% 6. Total input reactance, Xin
% 7. Length of the inner dipole, L'
% 8. Capacitance of each resonating capacitor, C
%
% The program also outputs (returns) the values of the paramaters input by the user:
% 1. Outer dipole radius, a
% 2. T-match wire radius, a'
% 3. Center-to-center wire separation, s
% 4. Outer dipole length, L ***OR*** the impedance of the antenna 
% 5. Characteristic impedance of the feed transmission line, Zc
% 6. Frequency at which to resonate, f
%-------------------------------------------------------------------------------------
% Written by Light Bryant-November 2002
%-------------------------------------------------------------------------------------
%
clear all;
close all;
format long;

%DEFINITION OF CONSTANTS
C = 0.577215665;
k = 2*pi;
eta = 376.73;

fprintf('Output device option \n\tOption (1): Screen\n\tOption (2): File \n');
ERR = 1;
while(ERR ~= 0)
   DEVICE = str2num(input('\nOutput device = ','s'));
   if(DEVICE == 1)
      ERR = 0;
   elseif(DEVICE == 2)
      FILNAM = input('Input the desired output filename: ','s');
      ERR = 0;
   else
      fprintf('\nOutput device number should be either 1 or 2\n');
   end
end

fprintf('Choose antenna option \n\tOption (3): Dipole of known length\n\tOption (4): Antenna of known impedance \n');
ERR = 1;
while(ERR ~= 0)
   BIFURCATE = str2num(input('\nAntenna option = ','s'));
   if(BIFURCATE == 3)
      ERR = 0;
   

%OBTAIN THE USER INPUTS
ERR = 1;
while ERR ~= 0
    a = str2num(input('\nThe wire radius of the dipole (in wavelengths) = ', 's'));
    a_prime = str2num(input('\nThe feed wire radius (in wavelengths) = ', 's'));
    s = str2num(input('\nThe center-to-center wire separation (in wavelengths) = ', 's'));
    if (a + a_prime) >= s
      fprintf('\n   *** ERROR: Please make sure the seperation is larger than the addition of the two radii!\n\n');
    else
      ERR =0;
    end
end
L = str2num(input('\nThe length of the longer (outer) dipole (in wavelengths) = ', 's'));
Xc = str2num(input('\nThe characteristic impedance (real) of the feed transmission line (in ohms) = ', 's'));
f = str2num(input('\nThe resonant (center) frequency (in megahertz) = ', 's'));


%CALCULATE THE CURRENT DIVISION FACTOR == ALPHA
u = a/a_prime;
v = s/a_prime;
alpha = (acosh((v^2-u^2+1)/(2*v)))/(acosh((v^2+u^2-1)/(2*v*u)));


%OBTAIN THE EFFECTIVE ("EQUIVALENT") RADIUS FOR THE GIVEN DIPOLE GEOMETRY
ln_ae = log(a_prime)+(1/((1+u)^2))*(u^2*log(u)+2*u*log(v));
ae = exp(ln_ae);


%DETERMINE THE INPUT IMPEDANCE OF THE (OUTER) DIPOLE IN THE ANTENNA (RADIATING) MODE
Rr = (eta)/(2*pi)*(C+log(k*L)-cosint(k*L)+0.5*sin(k*L)*(sinint(2*k*L)-2*sinint(k*L))+0.5*cos(k*L)*(C+log(0.5*k*L)+cosint(2*k*L)-2*cosint(k*L)));
Rina = Rr/((sin(0.5*k*L))^2);
Xm = (eta)/(4*pi)*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))-sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint((2*k*a^2)/L)));
Xina = Xm/((sin(0.5*k*L))^2);
Za = Rina + j*Xina;


%STEP UP ANTENNA MODE IMPEDANCE BY (1 + ALPHA)^2 TO PLACE IN SHUNT WITH TRANSMISSION LINE MODE; 
%DESIGNATE RESULT BY Z2 
Z2 = Za*(1+alpha)^2;


%DETERMINE THE CHARACTERISTIC IMPEDANCE AT THE INPUT TERMINALS FOR THE TRANSMISSION LINE MODE 
%(TWO-WIRE SHORTED TRANSMISSION LINE OF LENGTH L'/2 WITH GIVEN RADII AND SEPARATION DISTANCE)
Zo = (eta/k)*acosh((s^2-a^2-a_prime^2)/(2*a*a_prime));


%STEP-BY-STEP CALCULATION OF L'; STEPS INCLUDE NORMALIZING AND RENORMALIZING BY THE VALUE OF 
%THE IMPEDANCE, Zo, OF THE TRANSMISSION LINE FORMED BY THE DRIVEN ELEMENT AND THE T-MATCH 
%INNER DIPOLE.  THE ALGORITHM EFFECTIVELY WORKS BACKWARD BY APPLYING THE CONSTRAINT Rin = Zc 
%= Xc, WHICH THE T-MATCH IS DESIGNED TO SATISFY
z2 = Z2/Zo;
y2 = 1/z2;
g2 = real(y2); b2 = imag (y2);
bt_1 = -b2 + sqrt((Zo*g2 - Xc*g2^2) / Xc);
bt_2 = -b2 - sqrt((Zo*g2 - Xc*g2^2) / Xc);
L_1 = 2/k * acot(-2 * bt_1); L_2 = 2/k * acot(-2 * bt_2);
L_prime = L_2;	
Zt = j*Zo*tan(0.5*k*L_prime);
zt = Zt/Zo; 
yt = 1/(2*zt); gt = real(yt); bt = imag (yt);
%%%%%%formatting output
RT = real(Zt);
ZT = imag(Zt);
if ZT < 0
   sign_ZT = ' - ';
else
   sign_ZT = ' + ';
end
ZT = abs(ZT);
%%%%%%


%DETERMINE THE TOTAL INPUT IMPEDANCE, WHICH COMBINES THE ANTENNA (RADIATING) AND TRANSMISSION 
%(NONRADIATING) MODES
yin = y2 + yt;
zin = 1/yin;
rin = real(zin); xin = imag(zin);
Rin = Zo*rin; Xin = Zo*xin; Zin = Rin + j*Xin;
%%%%%%formatting output
RIN = real(Zin);
XIN = imag(Zin);
if XIN < 0
   sign_XIN = ' - ';
else
   sign_XIN = ' + ';
end
XIN = abs(XIN);
%%%%%%


%DETERMINE THE SERIES CAPACITANCE VALUE TO REMOVE THE REACTIVE IMPEDANCE AT THE INPUT TERMINALS
%FIRST CONVERT USER-INPUT FREQUENCY (in MHz) TO Hz
f_prime = f * 10^6; 
C = 1/(pi*f_prime*Xin);
Cin = C/2;


%FOR PURPOSES OF OUTPUT, IF THE TRANSMISSION LINE MODE IMPEDANCE IS EXTREMELY LARGE, ITS VALUE 
%IS DESIGNATED AS BEING INFINITE
if (abs(Zt) > 10^12)
   Zt = inf;
end


%DETERMINE WHETHER OUTPUT IS DIRECTED TO THE SCREEN OR TO A USER-NAMED FILE
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end


%DISPLAY INPUT PARAMETER VALUES AND CALCULATED OUTPUT RESULTS
fprintf(fid,'\n- - - - - - - - - - - - - - T-MATCH SUMMARY - - - - - - - - - - - - - -');
fprintf(fid,'\n\n-------------------------------INPUT DATA------------------------------');
fprintf(fid,'\n\nThe wire radius, a, of the dipole (in wavelengths) is %s\n',num2str(a));
fprintf(fid,'\nThe feed wire radius, aprime, (in wavelengths) is %s\n',num2str(a_prime));
fprintf(fid,'\nThe center-to-center wire separation, s, (in wavelengths) is %s\n',num2str(s));
fprintf(fid,'\nThe total length, L, of the dipole (in wavelengths) is %s\n',num2str(L));
fprintf(fid,'\nThe characteristic impedance, Zc, of the feed transmission line (in ohms) is %s\n',num2str(Xc));
fprintf(fid,'\nThe resonant (center) frequency (in megahertz) is %s\n',num2str(f));
fprintf(fid,'\n\n\n');
fprintf(fid,'---------------------------------OUTPUT--------------------------------');
fprintf(fid,'\n\nThe current division factor of the dipole pair is: alpha = %s\n',num2str(alpha));
fprintf(fid,'\nThe effective radius (in wavelengths) of the dipole pair is: ae = %s\n',num2str(ae));
fprintf(fid,'\nThe characteristic impedance (in ohms) of the two-dipole transmission line mode is: \nZo = %s\n', num2str(Zo));
fprintf(fid,'\nThe impedance (in ohms) at the input terminals for the two-dipole shorted transmission \nline mode is: Zt = %4.4f%sj%4.4f\n',RT,sign_ZT,ZT);
fprintf(fid,'\nThe total (radiating and nonradiating modes) two-dipole input impedance (in ohms) is: \nZin = %4.4f%sj%4.4f\n',RIN,sign_XIN,XIN);
fprintf(fid,'\nThe total dipole input resistance (in ohms) is: Rin = %s\n', num2str(Rin));
fprintf(fid,'\nThe total dipole input reactance (in ohms) is: Xin = %s\n', num2str(Xin));
fprintf(fid,'\nThe total input series capacitance (in farads) to resonate the input impedance is: \nCin = %s\n', num2str(Cin));
fprintf(fid,'\nThe individual capacitor series capacitance (in farads) to balance the line is: \nC = %s\n', num2str(C));
fprintf(fid,'\nThe inner dipole length (in wavelengths) is: Lprime = %s\n', num2str(L_prime));
fprintf(fid,'\n\n');

   elseif(BIFURCATE == 4)
      ERR = 0;
      
%OBTAIN THE USER INPUTS
a = str2num(input('\nThe wire radius of the antenna (in wavelengths) = ', 's'));
a_prime = str2num(input('\nThe T-match wire radius (in wavelengths) = ', 's'));
s = str2num(input('\nThe center-to-center wire separation (in wavelengths) = ', 's'));
ZAR = str2num(input('\nThe antenna input resistance (in ohms) = ', 's'));
ZAX = str2num(input('\nThe antenna input reactance (in ohms) = ', 's'));
Xc = str2num(input('\nThe characteristic impedance (real) of the feed transmission line (in ohms) = ', 's'));
f = str2num(input('\nThe resonant (center) frequency (in megahertz) = ', 's'));


%CALCULATE THE CURRENT DIVISION FACTOR == ALPHA
u = a/a_prime;
v = s/a_prime;
alpha = (acosh((v^2-u^2+1)/(2*v)))/(acosh((v^2+u^2-1)/(2*v*u)));


%OBTAIN THE EFFECTIVE ("EQUIVALENT") RADIUS FOR THE GIVEN GEOMETRY
ln_ae = log(a_prime)+(1/((1+u)^2))*(u^2*log(u)+2*u*log(v));
ae = exp(ln_ae);

%USE THE GIVEN ANTENNA IMPEDANCE
ZA = ZAR + j*ZAX;
Za = ZA;


%STEP UP ANTENNA MODE IMPEDANCE BY (1 + ALPHA)^2 TO PLACE IN SHUNT WITH TRANSMISSION LINE MODE; 
%DESIGNATE RESULT BY Z2 
Z2 = Za*(1+alpha)^2;
%%%%%%formatting output
ZAR = real(ZA);
ZAX = imag(ZA);
if ZAX < 0
   sign_ZAX = ' - ';
else
   sign_ZAX = ' + ';
end
ZAX = abs(ZAX);
%%%%%%


%DETERMINE THE CHARACTERISTIC IMPEDANCE AT THE INPUT TERMINALS FOR THE TRANSMISSION LINE MODE 
%(TWO-WIRE SHORTED TRANSMISSION LINE OF LENGTH L'/2 WITH GIVEN RADII AND SEPARATION DISTANCE)
Zo = (eta/k)*acosh((s^2-a^2-a_prime^2)/(2*a*a_prime));


%STEP-BY-STEP CALCULATION OF L'; STEPS INCLUDE NORMALIZING AND RENORMALIZING BY THE VALUE OF 
%THE IMPEDANCE, Zo, OF THE TRANSMISSION LINE FORMED BY THE DRIVEN ELEMENT AND THE T-MATCH 
%INNER DIPOLE.  THE ALGORITHM EFFECTIVELY WORKS BACKWARD BY APPLYING THE CONSTRAINT Rin = Zc 
%= Xc, WHICH THE T-MATCH IS DESIGNED TO SATISFY
z2 = Z2/Zo;
y2 = 1/z2; g2 = real(y2); b2 = imag (y2);
bt_1 = -b2 + sqrt((Zo*g2 - Xc*g2^2) / Xc);
bt_2 = -b2 - sqrt((Zo*g2 - Xc*g2^2) / Xc);
L_1 = 2/k * acot(-2 * bt_1); L_2 = 2/k * acot(-2 * bt_2);
L_prime = L_2;	
Zt = j*Zo*tan(0.5*k*L_prime);
zt = Zt/Zo; 
yt = 1/(2*zt); gt = real(yt); bt = imag (yt);
%%%%%%formatting output
RT = real(Zt);
ZT = imag(Zt);
if ZT < 0
   sign_ZT = ' - ';
else
   sign_ZT = ' + ';
end
ZT = abs(ZT);
%%%%%%


%DETERMINE THE TOTAL INPUT IMPEDANCE, WHICH COMBINES THE ANTENNA (RADIATING) AND TRANSMISSION 
%(NONRADIATING) MODES
yin = y2 + yt;
zin = 1/yin;
rin = real(zin); xin = imag(zin);
Rin = Zo*rin; Xin = Zo*xin; Zin = Rin + j*Xin;
%%%%%%formatting output
RIN = real(Zin);
XIN = imag(Zin);
if XIN < 0
   sign_XIN = ' - ';
else
   sign_XIN = ' + ';
end
XIN = abs(XIN);
%%%%%%


%DETERMINE THE SERIES CAPACITANCE VALUE TO REMOVE THE REACTIVE IMPEDANCE AT THE INPUT TERMINALS
%FIRST CONVERT USER-INPUT FREQUENCY (in MHz) TO Hz
f_prime = f * 10^6; 
C = 1/(pi*f_prime*Xin);
Cin = C/2;


%FOR PURPOSES OF OUTPUT, IF THE TRANSMISSION LINE MODE IMPEDANCE IS EXTREMELY LARGE, ITS VALUE 
%IS DESIGNATED AS BEING INFINITE
if (abs(Zt) > 10^12)
   Zt = inf;
end


%DETERMINE WHETHER OUTPUT IS DIRECTED TO THE SCREEN OR TO A USER-NAMED FILE
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end


%DISPLAY INPUT PARAMETER VALUES AND CALCULATED OUTPUT RESULTS
fprintf(fid,'\n- - - - - - - - - - - - - - T-MATCH SUMMARY - - - - - - - - - - - - - -');
fprintf(fid,'\n\n-------------------------------INPUT DATA------------------------------');
fprintf(fid,'\n\nThe wire radius, a, of the antenna (in wavelengths) is %s\n',num2str(a));
fprintf(fid,'\nThe T-match wire radius, aprime, (in wavelengths) is %s\n',num2str(a_prime));
fprintf(fid,'\nThe center-to-center wire separation, s, (in wavelengths) is %s\n',num2str(s));
fprintf(fid,'\nThe antenna input impedance, Za, (in ohms) is %4.4f%sj%4.4f\n',ZAR,sign_ZAX,ZAX);
fprintf(fid,'\nThe characteristic impedance, Zc, of the feed transmission line (in ohms) is %s\n',num2str(Xc));
fprintf(fid,'\nThe resonant (center) frequency (in megahertz) is %s\n',num2str(f));
fprintf(fid,'\n\n\n');
fprintf(fid,'---------------------------------OUTPUT--------------------------------');
fprintf(fid,'\n\nThe current division factor of the antenna pair is: alpha = %s\n',num2str(alpha));
fprintf(fid,'\nThe effective radius (in wavelengths) of the antenna pair is: ae = %s\n',num2str(ae));
fprintf(fid,'\nThe characteristic impedance (in ohms) of the two-element transmission line mode is: \nZo = %s\n', num2str(Zo));
fprintf(fid,'\nThe impedance (in ohms) at the input terminals for the two-dipole shorted transmission \nline mode is: Zt = %4.4f%sj%4.4f\n',RT,sign_ZT,ZT);
fprintf(fid,'\nThe total (radiating and nonradiating modes) two-dipole input impedance (in ohms) is: \nZin = %4.4f%sj%4.4f\n',RIN,sign_XIN,XIN);
fprintf(fid,'\nThe total input resistance (in ohms) is: Rin = %s\n', num2str(Rin));
fprintf(fid,'\nThe total input reactance (in ohms) is: Xin = %s\n', num2str(Xin));
fprintf(fid,'\nThe total input series capacitance (in farads) to resonate the input impedance is: \nCin = %s\n', num2str(Cin));
fprintf(fid,'\nThe individual capacitor series capacitance (in farads) to balance the line is: \nC = %s\n', num2str(C));
fprintf(fid,'\nThe inner T-match element length (in wavelengths) is: Lprime = %s\n', num2str(L_prime));
fprintf(fid,'\n\n');

   else
      fprintf('\nAntenna option choice should be either 3 or 4\n');
   end
end


if(DEVICE == 2)
   fclose(fid);
end

%---conclusion of program---------------------------------------------------------------------------------------
