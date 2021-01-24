%**************************************************************************************
% GAMMA.m
%*************************************************************************
% This is a MATLAB based program to analyze and design a Gamma_Match connection. 
%
% The gamma-match connection (which is a special case of the T-match, and thus a 
% particular case of the folded dipole) is employed as a broadband matching 
% technique to feed a dipole from a mismatched (and possibly unbalanced) 
% transmission line (such as a coaxial line).  The inner gamma rod (usually with 
% smaller radius) is connected to the transmission line at one end and is tapped at 
% one location to the outer dipole (or antenna), the tap point set a distance l'/2 
% from its center.  The gamma match and dipole are separated from each other by a 
% small distance, s.  A capacitor, whose purpose is to resonate the circuit and 
% balance the line, is placed at the feed point of the gamma rod.  (See Figure 9.21 
% from the book Antenna Theory by Balanis.)   
%
% This MATLAB program for a gamma-match connection is designed to calculate the values 
% of l'/2 and C, which are, respectively, the length of the gamma rod (fed by the
% transmission line) and the capacitance of the capacitor at the junction of the
% transmission line and the gamma rod.  
%
% The program assumes the geometry of the gamma rod-dipole is fixed (other than the 
% length of the gamma rod).  The user is prompted to select from one of two options: 
% EITHER a dipole antenna whose length is known OR an arbitrary antenna whose input 
% impedance is known.  For either case the user is prompted to input the values of the 
% two radii (whether the gamma rod plus dipole or the antenna with gamma-match connection) 
% and the separation between the two elements.  For the two-dipole case, the user inputs 
% the length of the outer dipole; for the antenna-T-match arrangement the user inputs the 
% impedance of the antenna (its real and imaginary parts separately).  
%
% In addition, the user inputs the characteristic impedance (which is assumed to be real-
% valued) of the feed transmission line; this value is denoted Zc.  The user must also 
% input the frequency (in megahertz) at which the dipole is to resonate.
%
%
%-------------------------------WARNING!!-------------------------------------------------
% THIS PROGRAM WILL FAIL TO YIELD ANY SATISFACTORY ANSWER IF THE TRANSMISSION LINE 
% CHARACTERISTIC IMPEDANCE, Zc, IS TOO LARGE.  IN PARTICULAR, IF Zc EXCEEDS Zo/g2, THE 
% CALCULATED OUTPUT LENGTH L'/2 WILL BE A COMPLEX NUMBER; THE VARIABLES Zo AND g2 ARE 
% CALCULATED WITHIN THE PROGRAM, AND THEIR VALUES DEPEND ON THE GEOMETRY SPECIFIED BY THE 
% USER (i.e., ON THE INPUT VALUES a, a', s, AND L).
%
% FURTHERMORE, IF Zc APPROACHES (FROM BELOW) TOO CLOSELY THE CRITICAL VALUE OF Zo/g2, THE 
% OUTPUT VALUE OF L'/2 MAY BE REAL, BUT IT WILL BE NEGATIVE.
%
% THE CENTER-TO-CENTER SEPARATION, s, OF THE GAMMA ROD AND THE OUTER DIPOLE MUST BE 
% SUFFICIENTLY LARGE TO SATISFY THE CONSTRAINT s^2 > (a^2 + (a')^2).  THIS CONSTRAINT SHOULD 
% NEVER POSE A PROBLEM FOR REALISTIC GEOMETRIES, BUT MUST BE HEEDED TO PREVENT THE USER FROM 
% ENTERING TOO SMALL A VALUE OF s.
%-----------------------------------------------------------------------------------------
%
%
% The program computes and outputs the following values:
% the effective radius of the gamma rod-dipole configuration, ae
% the dipole (antenna) input impedance, Za
% the characteristic impedance of the element pair, Zo
% the transmission line mode impedance, Zt
% the total input resistance, Rin
% the total input reactance, Xin
% the length of the gamma rod, L'/2 = Lgamma
% the capacitance of the resonating capacitor, C
%
% The program also outputs (returns) the values of the paramaters input by the user:
% the outer dipole (antenna) radius, a
% the gamma-match wire radius, a'
% the center-to-center wire separation, s
% the outer dipole length, L ***OR*** the impedance of the antenna 
% the characteristic impedance of the feed transmission line, Zc
% the frequency at which to resonate, f
%----------------------------------------------------------------------------------
% Written by Light Bryant-November 2002 
%----------------------------------------------------------------------------------
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
a = str2num(input('\nThe wire radius of the dipole (in wavelengths) = ', 's'));
a_prime = str2num(input('\nThe feed wire radius (in wavelengths) = ', 's'));
s = str2num(input('\nThe center-to-center wire separation (in wavelengths) = ', 's'));
L = str2num(input('\nThe length of the dipole (in wavelengths) = ', 's'));
Xc = str2num(input('\nThe characteristic impedance (real) of the feed transmission line (in ohms) = ', 's'));
f = str2num(input('\nThe resonant (center) frequency (in megahertz) = ', 's'));


%CALCULATE THE CURRENT DIVISION FACTOR == ALPHA
u = a/a_prime;
v = s/a_prime;
alpha = (acosh((v^2-u^2+1)/(2*v)))/(acosh((v^2+u^2-1)/(2*v*u)));


%OBTAIN THE EFFECTIVE ("EQUIVALENT") RADIUS FOR THE GIVEN TWO-ELEMENT GEOMETRY
ln_ae = log(a_prime)+(1/((1+u)^2))*(u^2*log(u)+2*u*log(v));
ae = exp(ln_ae);


%DETERMINE THE INPUT IMPEDANCE OF THE DIPOLE IN THE ANTENNA (RADIATING) MODE
Rr = (eta)/(2*pi)*(C+log(k*L)-cosint(k*L)+0.5*sin(k*L)*(sinint(2*k*L)-2*sinint(k*L))+0.5*cos(k*L)*(C+log(0.5*k*L)+cosint(2*k*L)-2*cosint(k*L)));
Rina = Rr/((sin(0.5*k*L))^2);
Xm = (eta)/(4*pi)*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))-sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint((2*k*a^2)/L)));
Xina = Xm/((sin(0.5*k*L))^2);
Za = Rina + j*Xina;


%STEP UP ANTENNA MODE IMPEDANCE BY (1 + ALPHA)^2 TO PLACE IN SHUNT WITH TRANSMISSION LINE MODE; 
%DESIGNATE RESULT BY Z2 
Z2 = Za/2*(1+alpha)^2;


%DETERMINE THE CHARACTERISTIC IMPEDANCE AT THE INPUT TERMINALS FOR THE TRANSMISSION LINE MODE 
%(TWO-WIRE SHORTED TRANSMISSION LINE OF LENGTH L'/2 WITH GIVEN RADII AND SEPARATION DISTANCE)
Zo = (eta/k)*acosh((s^2-a^2-a_prime^2)/(2*a*a_prime));


%STEP-BY-STEP CALCULATION OF L'/2; STEPS INCLUDE NORMALIZING AND RENORMALIZING BY THE VALUE OF 
%THE IMPEDANCE, Zo, OF THE TRANSMISSION LINE FORMED BY THE DRIVEN ELEMENT AND THE T-MATCH 
%INNER DIPOLE.  THE ALGORITHM EFFECTIVELY WORKS BACKWARD BY APPLYING THE CONSTRAINT Rin = Zc 
%= Xc, WHICH THE T-MATCH IS DESIGNED TO SATISFY
z2 = Z2/Zo;
y2 = 1/z2; g2 = real(y2); b2 = imag (y2);
bt_1 = -b2 + sqrt((Zo*g2 - Xc*g2^2) / Xc);
bt_2 = -b2 - sqrt((Zo*g2 - Xc*g2^2) / Xc);
L_1 = 2/k * acot(-1 * bt_1); L_2 = 2/k * acot(-1 * bt_2);
L_prime = L_2;
L_gamma = L_prime/2;
Zg = j*Zo*tan(0.5*k*L_prime);
zg = Zg/Zo; 
yg = 1/(1*zg); gg = real(yg); bg = imag (yg);
%%%%%%formatting output
RG = real(Zg);
XG = imag(Zg);
if XG < 0
   sign_XG = ' - ';
else
   sign_XG = ' + ';
end
XG = abs(XG);
%%%%%%


%DETERMINE THE TOTAL INPUT IMPEDANCE, WHICH COMBINES THE ANTENNA (RADIATING) AND TRANSMISSION 
%(NONRADIATING) MODES
yin = y2 + yg;
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
C = 1/(2*pi*f_prime*Xin);


%FOR PURPOSES OF OUTPUT, IF THE TRANSMISSION LINE MODE IMPEDANCE IS EXTREMELY LARGE, ITS VALUE 
%IS DESIGNATED AS BEING INFINITE
if (abs(Zg) > 10^12)
   Zg = inf;
end


%DETERMINE WHETHER OUTPUT IS DIRECTED TO THE SCREEN OR TO A USER-NAMED FILE
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end


%DISPLAY INPUT PARAMETER VALUES AND CALCULATED OUTPUT RESULTS
fprintf(fid,'\n- - - - - - - - - - - - - GAMMA-MATCH SUMMARY - - - - - - - - - - - - -');
fprintf(fid,'\n\n-------------------------------INPUT DATA------------------------------');
fprintf(fid,'\n\nThe wire radius, a, of the dipole (in wavelengths) is %s\n',num2str(a));
fprintf(fid,'\nThe gamma rod radius, aprime, (in wavelengths) is %s\n',num2str(a_prime));
fprintf(fid,'\nThe center-to-center wire separation, s, (in wavelengths) is %s\n',num2str(s));
fprintf(fid,'\nThe total length, L, of the dipole (in wavelengths) is %s\n',num2str(L));
fprintf(fid,'\nThe characteristic impedance, Zc, of the feed transmission line (in ohms) is %s\n',num2str(Xc));
fprintf(fid,'\nThe resonant (center) frequency (in megahertz) is %s\n',num2str(f));
fprintf(fid,'\n\n\n');
fprintf(fid,'---------------------------------OUTPUT--------------------------------');
fprintf(fid,'\n\nThe current division factor is: alpha = %s\n',num2str(alpha));
fprintf(fid,'\nThe effective radius (in wavelengths) is: ae = %s\n',num2str(ae));
fprintf(fid,'\nThe characteristic impedance (in ohms) of the gamma rod-dipole transmission line mode is: \nZo = %s\n', num2str(Zo));
fprintf(fid,'\nThe impedance (in ohms) at the input terminals for the two-dipole shorted transmission \nline mode is: Zg = %4.4f%sj%4.4f\n',RG,sign_XG,XG);
fprintf(fid,'\nThe total (radiating and nonradiating modes) two-dipole input impedance (in ohms) is: \nZin = %4.4f%sj%4.4f\n',RIN,sign_XIN,XIN);
fprintf(fid,'\nThe total input resistance (in ohms) is: Rin = %s\n', num2str(Rin));
fprintf(fid,'\nThe total input reactance (in ohms) is: Xin = %s\n', num2str(Xin));
fprintf(fid,'\nThe series capacitance (in farads) to resonate the input impedance and to balance the \nline is: C = %s\n', num2str(C));
fprintf(fid,'\nThe gamma rod length (in wavelengths) is: Lgamma = %s\n', num2str(L_gamma));
fprintf(fid,'\n\n');

   elseif(BIFURCATE == 4)
      ERR = 0;
      
      
%OBTAIN THE USER INPUTS
a = str2num(input('\nThe wire radius of the dipole (in wavelengths) = ', 's'));
a_prime = str2num(input('\nThe feed wire radius (in wavelengths) = ', 's'));
s = str2num(input('\nThe center-to-center wire separation (in wavelengths) = ', 's'));
ZAR = str2num(input('\nThe antenna input resistance (in ohms) = ', 's'));
ZAX = str2num(input('\nThe antenna input reactance (in ohms) = ', 's'));
Xc = str2num(input('\nThe characteristic impedance (real) of the feed transmission line (in ohms) = ', 's'));
f = str2num(input('\nThe resonant (center) frequency (in megahertz) = ', 's'));


%CALCULATE THE CURRENT DIVISION FACTOR == ALPHA
u = a/a_prime;
v = s/a_prime;
alpha = (acosh((v^2-u^2+1)/(2*v)))/(acosh((v^2+u^2-1)/(2*v*u)));


%OBTAIN THE EFFECTIVE ("EQUIVALENT") RADIUS FOR THE GIVEN TWO-ELEMENT GEOMETRY
ln_ae = log(a_prime)+(1/((1+u)^2))*(u^2*log(u)+2*u*log(v));
ae = exp(ln_ae);
      
      
%USE THE GIVEN ANTENNA IMPEDANCE
ZA = ZAR + ZAX*i;
Za = ZA;


%STEP UP ANTENNA MODE IMPEDANCE BY (1 + ALPHA)^2 TO PLACE IN SHUNT WITH TRANSMISSION LINE MODE; 
%DESIGNATE RESULT BY Z2 
Z2 = Za/2*(1+alpha)^2;
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


%STEP-BY-STEP CALCULATION OF L'/2; STEPS INCLUDE NORMALIZING AND RENORMALIZING BY THE VALUE OF 
%THE IMPEDANCE, Zo, OF THE TRANSMISSION LINE FORMED BY THE DRIVEN ELEMENT AND THE T-MATCH 
%INNER DIPOLE.  THE ALGORITHM EFFECTIVELY WORKS BACKWARD BY APPLYING THE CONSTRAINT Rin = Zc 
%= Xc, WHICH THE T-MATCH IS DESIGNED TO SATISFY
z2 = Z2/Zo;
y2 = 1/z2; g2 = real(y2); b2 = imag (y2);
bt_1 = -b2 + sqrt((Zo*g2 - Xc*g2^2) / Xc);
bt_2 = -b2 - sqrt((Zo*g2 - Xc*g2^2) / Xc);
L_1 = 2/k * acot(-1 * bt_1); L_2 = 2/k * acot(-1 * bt_2);
L_prime = L_2;
L_gamma = L_prime/2;
Zg = j*Zo*tan(0.5*k*L_prime);
zg = Zg/Zo; 
yg = 1/(1*zg); gg = real(yg); bg = imag (yg);
%%%%%%formatting output
RG = real(Zg);
XG = imag(Zg);
if XG < 0
   sign_XG = ' - ';
else
   sign_XG = ' + ';
end
XG = abs(XG);
%%%%%%


%DETERMINE THE TOTAL INPUT IMPEDANCE, WHICH COMBINES THE ANTENNA (RADIATING) AND TRANSMISSION 
%(NONRADIATING) MODES
yin = y2 + yg;
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
C = 1/(2*pi*f_prime*Xin);


%FOR PURPOSES OF OUTPUT, IF THE TRANSMISSION LINE MODE IMPEDANCE IS EXTREMELY LARGE, ITS VALUE 
%IS DESIGNATED AS BEING INFINITE
if (abs(Zg) > 10^12)
   Zg = inf;
end


%DETERMINE WHETHER OUTPUT IS DIRECTED TO THE SCREEN OR TO A USER-NAMED FILE
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end


%DISPLAY INPUT PARAMETER VALUES AND CALCULATED OUTPUT RESULTS
fprintf(fid,'\n- - - - - - - - - - - - - GAMMA-MATCH SUMMARY - - - - - - - - - - - - -');
fprintf(fid,'\n\n-------------------------------INPUT DATA------------------------------');
fprintf(fid,'\n\nThe wire radius, a, of the antenna (in wavelengths) is %s\n',num2str(a));
fprintf(fid,'\nThe gamma rod radius, aprime, (in wavelengths) is %s\n',num2str(a_prime));
fprintf(fid,'\nThe center-to-center wire separation, s, (in wavelengths) is %s\n',num2str(s));
fprintf(fid,'\nThe antenna input impedance, Za, (in ohms) is %4.4f%sj%4.4f\n',ZAR,sign_ZAX,ZAX);
fprintf(fid,'\nThe characteristic impedance, Zc, of the feed transmission line (in ohms) is %s\n',num2str(Xc));
fprintf(fid,'\nThe resonant (center) frequency (in megahertz) is %s\n',num2str(f));
fprintf(fid,'\n\n\n');
fprintf(fid,'---------------------------------OUTPUT--------------------------------');
fprintf(fid,'\n\nThe current division factor is: alpha = %s\n',num2str(alpha));
fprintf(fid,'\nThe effective radius (in wavelengths) is: ae = %s\n',num2str(ae));
fprintf(fid,'\nThe characteristic impedance (in ohms) of the gamma rod-dipole transmission line mode is: \nZo = %s\n', num2str(Zo));
fprintf(fid,'\nThe impedance (in ohms) at the input terminals for the two-dipole shorted transmission \nline mode is: Zg = %4.4f%sj%4.4f\n',RG,sign_XG,XG);
fprintf(fid,'\nThe total (radiating and nonradiating modes) two-dipole input impedance (in ohms) is: \nZin = %4.4f%sj%4.4f\n',RIN,sign_XIN,XIN);
fprintf(fid,'\nThe total input resistance (in ohms) is: Rin = %s\n', num2str(Rin));
fprintf(fid,'\nThe total input reactance (in ohms) is: Xin = %s\n', num2str(Xin));
fprintf(fid,'\nThe series capacitance (in farads) to resonate the input impedance and to balance the \nline is: C = %s\n', num2str(C));
fprintf(fid,'\nThe gamma rod length (in wavelengths) is: Lgamma = %s\n', num2str(L_gamma));
fprintf(fid,'\n\n');

   else
      fprintf('\nAntenna option choice should be either 3 or 4\n');
   end
end


if(DEVICE == 2)
   fclose(fid);
end

%---conclusion of program---------------------------------------------------------------------------------------
