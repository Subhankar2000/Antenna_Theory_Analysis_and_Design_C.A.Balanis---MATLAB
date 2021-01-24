%***************************************************************************************
% FOLDED_DIPOLE.m
%************************************************************************
% This is a MATLAB based program written to aid in the design of folded dipole antennas.
%
% The folded dipole is used in applications where the characteristic
% impedance, Zc, of the input transmission line is larger than 50 or 75
% ohms.  In these situations, the single element dipole would not provide
% good matching characteristics with the transmission line.  In order to
% provide better matching characteristics, a folded wire, which forms a very
% thin (s<<lambda) rectangular loop, is used.  To generalize the design, the
% radii of the two wires of the folded dipole are not necessarily the same.
% The radius of the dipole is represented by a and that of the feed wire
% element is a'. The folded dipole serves as a step-up impedance
% transformer.  This is very useful when used in conjuction with high
% impedance balanced transmission lines such as "twin-lead" transmission
% lines.  Two capacitors, C, of equal value are placed in series at the feed
% point of the dipole to allow the circuit to resonate and to keep the line
% balanced.  In the case that the total input impedance is capacitive, a
% single inductor, Lin, is placed in parallel at the feed point of the
% dipole to allow the circuit to resonate and to keep the line balanced.
%
% To analyze this antenna, the user must enter the following quantities:
%       1.  Dipole radius, a, (in wavelengths)
%       2.  Feed wire dipole radius, a', (in wavelengths)
%       3.  Center-to-center wire separation, s, (in wavelengths)
%       If impedance is calculated for dipole of known length
%           4.  Dipole length, L, (in wavelengths)
%       If antenna impedance is specified
%           5.  The antenna input resistance, Rd, (in ohms)
%           6.  The antenna input reactance, Xd, (in ohms)
%       7.  Frequency for which to resonate, f, (in MHz)
%       8.  The real characteristic impedance of the line feeding the antenna, Zc, (in ohms)
%
% The program computes the following quantities:
%       1.  The wire current division factor, alpha
%       2.  Effective radius, ae, (in wavelengths)
%       3.  Input impedance of the single element dipole, Zd, (in ohms)
%       4.  Characteristic impedance for the transmission line mode, Zo, 
%           (in ohms)
%       5.  Transmission line mode impedance, Zt, (in ohms)
%       6.  Total input resistance, Rin, (in ohms)
%       7.  Total input reactance, Xin, (in ohms)
%       If the total input reactance is inductive
%           8.  Total input capacitance, Cin, to resonate the dipole 
%           (in farads)
%           9.  Capacitance of each resonating capacitor, C, to balance the 
%           line(in farads)
%       If the total input reactance is capacitive
%           10. Total input inductance, Lin, to resonate the dipole (in
%           henries)
%       11. Reflection coefficient after the total input impedance is
%           resonated
%       12. VSWR after the total input impedance is resonated
%-----------------------------------------------------------------------------------------------
%  Written by Jason Latimer, November 2002
%-----------------------------------------------------------------------------------------------
clear all;
close all;
format long;

%Definition of constants
C = 0.5772;
k = 2*pi;

%Flag variable
cap_used = 0;

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

%Prompt user to determine if antenna impedance will be calculated or given
fprintf('\n\nAntenna impedance option \n\tOption (3): Calculate impedance for dipole of known length\n\tOption (4): Specify total folded dipole impedance \n');
ERR = 1;
while(ERR ~= 0)
   IMP_OPT = str2num(input('\nAntenna impedance option = ','s'));
   if(IMP_OPT == 3)
      ERR = 0;
   elseif(IMP_OPT == 4)
      ERR = 0;
   else
      fprintf('\nAntenna impedance option should be either 3 or 4\n');
   end
end

%Prompt user to determine if shunt or series matching will be used
fprintf('\n\nAntenna matching option \n\tOption (5): Series matching\n\tOption (6): Shunt matching\n');
ERR = 1;
while(ERR ~= 0)
   MAT_OPT = str2num(input('\nAntenna matching option = ','s'));
   if(MAT_OPT == 5)
      ERR = 0;
   elseif(MAT_OPT == 6)
      ERR = 0;
   else
      fprintf('\nAntenna matching option should be either 5 or 6\n');
   end
end

%Obtain user inputs
if (IMP_OPT == 3)
    a = str2num(input('\nThe radius, a, of the dipole (in wavelengths) = ', 's'));
    a_prime = str2num(input('\nThe feed wire radius, aprime, (in wavelengths) = ', 's'));
    s = str2num(input('\nThe center-to-center wire separation, s, (in wavelengths) = ', 's'));
    L = str2num(input('\nThe total length, L, of the dipole (in wavelengths) = ', 's'));
end
if (IMP_OPT == 4)
    Rin_tot = str2num(input('\nThe folded dipole input resistance, Rin, (in ohms) = ', 's'));
    Xin_tot = str2num(input('\nThe folded dipole input reactance, Xin, (in ohms) = ', 's'));
    Zin_tot = Rin_tot + j*Xin_tot;
end
f = str2num(input('\nThe center frequency, f, (in MHz) = ', 's'));
Zc = str2num(input('\nThe real characteristic impedance of the line feeding the antenna, Zc, (in ohms) = ', 's'));

if (IMP_OPT == 3)
    %Calculate the current division factor == alpha
    u = a/a_prime;
    v = s/a_prime;
    alpha = (acosh((v^2-u^2+1)/(2*v)))/(acosh((v^2+u^2-1)/(2*v*u)));

    %Obtain the effective radius for the given dipole geometry
    ln_ae = log(a_prime)+(1/((1+u)^2))*(u^2*log(u)+2*u*log(v));
    ae = exp(ln_ae);

    %Determine the input impedance of the (outer) dipole in the antenna
    %(radiating) mode
    %note: This method is only used if option three is chosen
    Rr = (120*pi)/(2*pi)*(C+log(k*L)-cosint(k*L)+0.5*sin(k*L)*(sinint(2*k*L)-2*sinint(k*L))+0.5*cos(k*L)*(C+log(0.5*k*L)+cosint(2*k*L)-2*cosint(k*L)));
    Rin = Rr/((sin(0.5*k*L))^2);
    Xm = (120*pi)/(4*pi)*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))-sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint((2*k*a^2)/L)));
    Xin = Xm/((sin(0.5*k*L))^2);
    Zin = Rin + j*Xin;
    Za = Zin;
    Ra = real(Za);
    Xa = imag(Za);
    
    %Determine the characteristic impedance at the input terminals for the
    %transmission line mode
    Zo = 60*acosh((s^2-a^2-a_prime^2)/(2*a*a_prime));
    Zt = j*Zo*tan(0.5*k*L);
    Rt = real(Zt);
    Xt = imag(Zt);

    %Determine the total input impedance, which combines the antenna
    %(radiating) and transmission line modes
    Zin_tot = 2*Zt*(Za*(1+alpha)^2)/(2*Zt+Za*(1+alpha)^2);
    Rin_tot = real(Zin_tot);
    Xin_tot = imag(Zin_tot);
end

%Portion of coded dedicated to series matching
if (MAT_OPT == 5)
    X_LE = - Xin_tot;
    if (X_LE < 0)
        %Determine the series capacitance to resonate the circuit
        C = 1/(pi*(f*(10^6))*Xin_tot);
        Cin = C/2;
        cap_used = 1;
    else
        %Determine the series inductance to resonate the circuit (need to
        %verify this section)
        Lin = X_LE/(2*pi*f*10^6);
        L = Lin/2;
    end
    Rin_tot_new = Rin_tot;
end

%Portion of code dedicated to shunt matching
if (MAT_OPT == 6)
    Yin_tot = 1/Zin_tot;
    Bin_tot = imag(Yin_tot);
    B_L = - j*Bin_tot;
    X_LE = 1/B_L;
    if (imag(X_LE) < 0)
        %Determine the shunt capacitance to resonate the circuit
        Cin = 1/(j*X_LE*2*pi*f*10^6);
        cap_used = 1;
    else
        %Determine the shunt inductance to resonate the circuit
        Lin = X_LE/(j*2*pi*f*10^6);
    end
    Rin_tot_new = 1/(real(Yin_tot));
end

%Determine the VSWR and the reflection coefficient
if (MAT_OPT == 6)
    gamma = (Rin_tot_new - Zc)/(Rin_tot_new + Zc);
else
    gamma = (Rin_tot - Zc)/(Rin_tot + Zc);
end
VSWR = (1 + abs(gamma))/(1 - abs(gamma));

%Determine if the output is directed to the screen or a user-named file
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end

%Display input parameter values and calculated output results
fprintf(fid,'\n\n- - - - - - - - - - - - - Folded Dipole - - - - - - - - - - - - - - - -');
fprintf(fid,'\n------------------------Summary of input data--------------------------');
if (IMP_OPT == 3)
    fprintf(fid,'\nThe radius, a, of the dipole (in wavelengths) is %s\n',num2str(a));
    fprintf(fid,'\nThe feed wire radius, aprime, (in wavelengths) is %s\n',num2str(a_prime));
    fprintf(fid,'\nThe center-to-center wire separation, s, (in wavelengths) is %s\n',num2str(s))
    fprintf(fid,'\nThe total length, L, of the dipole (in wavelengths) is %s\n',num2str(L));
elseif (IMP_OPT == 4)
    fprintf(fid,'\nThe folded dipole input resistance, Rin, (in ohms) is %s\n',num2str(Rin_tot));
    fprintf(fid,'\nThe folded dipole input reactance, Xin, (in ohms) is %s\n',num2str(Xin_tot));
end
fprintf(fid,'\nThe center frequency, f, (in MHz) is %s\n',num2str(f));
fprintf(fid,'\nThe real characteristic impedance of the line feeding the antenna, Zc, (in ohms) is %s\n', num2str(Zc));
fprintf(fid,'\n\n');
fprintf(fid,'----------------------------------Output----------------------------------');
if (IMP_OPT == 3)
    fprintf(fid,'\nThe current division ratio: alpha = %s\n',num2str(alpha));
    fprintf(fid,'\nThe effective radius (in wavelengths) of the dipole pair: ae = %s\n',num2str(ae));
    fprintf(fid,'\nThe dipole input resistance (in ohms): Rd = %s\n', num2str(Ra));
    fprintf(fid,'\nThe dipole input reactance (in ohms): Xd = %s\n', num2str(Xa));
    fprintf(fid,'\nThe characteristic impedance (in ohms) for the transmission line mode: Zo = %s\n', num2str(Zo));
    if abs(Zt) < (10^12)
        fprintf(fid,'\nThe resistance (in ohms) at the input terminals for the transmission line mode: Rt = %s\n', num2str(Rt));
        fprintf(fid,'\nThe reactance (in ohms) at the input terminals for the transmission line mode: Xt = %s\n', num2str(Xt));
    else
        fprintf(fid,'\nThe impedance (in ohms) at the input terminals for the transmission line mode: Zt = infinity\n');
    end
end
fprintf(fid,'\nThe total input resistance (in ohms): Rin = %s\n', num2str(Rin_tot));
fprintf(fid,'\nThe total input reactance (in ohms): Xin = %s\n', num2str(Xin_tot));
if (MAT_OPT == 5)
    if (cap_used == 1)
        fprintf(fid,'\nThe total input series capacitance (in farads) to resonate the input impedance: Cin = %s\n', num2str(Cin));
        fprintf(fid,'\nThe input series capacitance (in farads) to balance the line: C = %s\n', num2str(C));
    else
        fprintf(fid,'\nThe total input series inductance (in henries) to resonate the input impedance: Lin = %s\n', num2str(Lin));
        fprintf(fid,'\nThe input series inductance (in henries) to balance the line: L = %s\n', num2str(L));
    end
end
if (MAT_OPT == 6)
    if (cap_used == 1)
        fprintf(fid,'\nThe total input parallel capacitance (in farads) to resonate the input impedance: Cin = %s\n', num2str(Cin));
    else
        fprintf(fid,'\nThe total input parallel inductance (in henries) to resonate the input impedance: Lin = %s\n', num2str(Lin));
    end
end
fprintf(fid,'\nThe new input resistance after resonance is: Rin(total) = %s\n', num2str(Rin_tot_new));
fprintf(fid,'\nThe reflection coefficient after the impedance is resonated: gamma = %s\n', num2str(gamma));
fprintf(fid,'\nThe VSWR after the impedance is resonated is: VSWR = %s\n', num2str(VSWR));
fprintf(fid,'\n\n');

if(DEVICE == 2)
   fclose(fid);
end

%---End of program-------------------------------------------
