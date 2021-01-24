%***********************************************************************
%   IMPEDANCE.m
%********************************************************************** 
%
%     SELF AND MUTUAL IMPEDANCES
%******************************************************************
%     THIS IS A MATLAB BASED PROGRAM THAT COMPUTES:
%
%        I.  SELF IMPEDANCE OF ANY LENGTH DIPOLE
%        II. MUTUAL IMPEDANCE BETWEEN TWO IDENTICAL LINEAR DIPOLES
%
%     BASED ON THE INDUCED EMF METHOD AND THE IDEAL CURRENT DISTRIBUTION
%     OF EQUATION (4-56).
%
%        I.  SELF IMPEDANCE (INPUT IMPEDANCE)
%
%            1.  BASED ON CURRENT AT THE INPUT
%                  Zin = Rin + jXin
%                  Rin = INPUT RESISTANCE [EQU. (8-60a) & (8-61a)]
%                  Xin = INPUT REACTANCE [EQU. (8-60b) & (8-61b)]
%            2.  BASED ON CURRENT MAXIMUM
%
%                  Zinm = Rinm + jXinm
%                  Rinm = SELF RESISTANCE [EQU. (8-60a)]
%                  Xinm = SELF REACTANCE [EQU. (8-60b)]
%
%        II.  MUTUAL IMPEDANCE 
%
%            1.  BASED ON CURRENT AT THE INPUT
%                  Z21i = R21i + jX21i [EQU. (8-68)]
%                  R21i = INPUT MUTUAL RESISTANCE 
%                  Xin = INPUT MUTUAL REACTANCE 
%            2. BASED ON CURRENT MAXIMUM
%                  Z21m = R21m + jX21m [EQU. (8-70)]
%                  R21m = SELF MUTUAL RESISTANCE
%                  X21m = SELF MUTUAL REACTANCE
%
%     THE DIPOLES FOR MUTUAL IMPEDANCE COMPUTATIONS MUST BE IDENTICAL
%     WITH LENGTH OF ODD MULTIPLES OF HALF WAVELENGTH. 
%
%        OPTION I.  SELF IMPEDANCE (INPUT IMPEDANCE)       
%     
%          ** ARRAY INPUT PARAMETERS
%
%          1.  LENGTH OF DIPOLE (IN WAVELENGTHS)
%          2.  RADIUS OF DIPOLE (IN WAVELENGTHS)
%
%       OPTION II.  MUTUAL IMPEDANCE
%
%         CHOICE A:  SIDE-BY-SIDE [Fig. 8.19(a)]
%
%           ** ARRAY INPUT PARAMETERS:
%
%           1.  LENGTH OF THE DIPOLES (IN WAVELENGTHS)
%           2.  HORIZONTAL DISPLACEMENT OF DIPOLES (IN WAVELENGTHS)
%
%         CHOICE B: COLLINEAR [Fig. 8.19(b)]
%
%          ** ARRAY INPUT PARAMETERS:
%
%          1.  LENGTH OF THE DIPOLES (IN WAVELENGTHS)
%          2.  VERTICAL DISPLACEMENT OF DIPOLES (IN WAVELENGTHS)
%
%         CHOICE C. PARALLEL-IN-ECHELON [Fig. 8.19(c)]
%
%          ** ARRAY INPUT PARAMETERS:
%
%          1.  LENGTH OF THE DIPOLES (IN WAVELENGTHS)
%          2.  VERTICAL DISPLACEMENT OF DIPOLES (IN WAVELENGTHS)
%          3.  HORIZONTAL DISPLACEMENT OF DIPOLES (IN WAVELENGTHS)
%
%**   NOTE: ALL THE INPUT PARAMETERS ARE IN WAVELENGTHS.
%******************************************************************
%     Written by: Panayiotis Anastasiou, Arizona State University
%******************************************************************
clear all;
close all;
k = 2*pi;
ETA = 120*pi;
C=0.5772156649;

% Choose Output Device --->  Screen or Output File 

device_option = sprintf('OUTPUT DEVICE OPTION\n');
device_option = sprintf('%s   OPTION (1): SCREEN\n',device_option);
device_option = sprintf('%s   OPTION (2): OUTPUT FILE\n\n',device_option);
device_option = sprintf('%sOUTPUT DEVICE =  ',device_option);

output_device = input(device_option);
if output_device == 2
   filename = sprintf('\nINPUT THE DESIRED OUTPUT FILENAME (in single quotes) = ');
   filename = input(filename);
end


    impedance_option = sprintf('\nSELF OR MUTUAL IMPEDANCE OPTION\n');
    impedance_option = sprintf('%s   OPTION (1): SELF IMPEDANCE\n',impedance_option);
    impedance_option = sprintf('%s   OPTION (2): MUTUAL IMPEDANCE\n\n',impedance_option);
    impedance_option = sprintf('%sOPTION NUMBER = ',impedance_option);

    option_number = input(impedance_option);

        switch option_number
        
            case 1
        
            % For self Impedance

            selfimpedance_length = sprintf('\n   LENGTH OF THE DIPOLE (IN WAVELENGTHS) = ');
            length_dipole = input(selfimpedance_length);
     
            selfimpedance_radius = sprintf('   RADIUS OF THE DIPOLE (IN WAVELENGTHS) = ');
            radius_dipole = input(selfimpedance_radius);

            l = length_dipole;
            a = radius_dipole;
           
        % Self impedance calculations of a dipole

            % Based on current maximum

            Rinm=(ETA/(2*pi))*(C+log(k*l)-cosint(k*l)+0.5*sin(k*l)*(sinint(2*k*l)-2*...
                 sinint(k*l))+0.5*cos(k*l)*(C+log(k*l/2)+cosint(2*k*l)-2*cosint(k*l)));

            Xinm=(ETA/(4*pi))*(2*sinint(k*l)+cos(k*l)*(2*sinint(k*l)-sinint(2*k*l))...
                 -sin(k*l)*(2*cosint(k*l)-cosint(2*k*l)-cosint((2*k*a^2)/l)));
 
            % Based on current at the input
           
            Rin = Rinm/(sin(k*l/2))^2;       
            Xin = Xinm/(sin(k*l/2))^2; 
    
        if output_device == 2
            diary(filename);
        end
     

            disp(' '),
            disp('   ************PROGRAM OUTPUT************'),
            disp('   SELF IMPEDANCE CALCULATION OF A DIPOLE'),
            disp(' '),
            disp([selfimpedance_length num2str(length_dipole,'%12.6f')]), 
            disp([selfimpedance_radius num2str(radius_dipole,'%12.6f')]), 
            disp(' '),
           
            if Rin > 10^15
                disp('   Rin  =     INFINITY')
            else            
            disp(['   Rin  =     '  num2str(Rin,'%12.6f')]),
            end
            
            if Rin > 10^15
                disp('   Xin  =     INFINITY')
            else
            disp(['   Xin  =     '  num2str(Xin,'%12.6f')]),
            end
            
            disp(['   Rinm =     '  num2str(Rinm,'%12.6f')]),
            disp(['   Xinm =     '  num2str(Xinm,'%12.6f')]), 
	 
          diary off;
           
            case 2

            % For mutual impedance

            dipoleconf_option = sprintf('\nDIPOLE CONFIGURATION OPTION\n');
            dipoleconf_option = sprintf('%s   OPTION (A): SIDE-BY-SIDE\n',dipoleconf_option);
            dipoleconf_option = sprintf('%s   OPTION (B): COLINEAR\n',dipoleconf_option);
            dipoleconf_option = sprintf('%s   OPTION (C): PARALLEL-IN-ECHELLON\n',dipoleconf_option);
            dipoleconf_option = sprintf('%s\nCONFIGURATION OPTION: A, B, or C',dipoleconf_option);
            dipoleconf_option = sprintf('%s (in single quotes) = ',dipoleconf_option);

            configuration_option = input(dipoleconf_option,'s');   

                switch configuration_option
                case {'''A''', '''a'''}

                    % If user inserts option A


                    optionA_length = sprintf('\n   LENGTH OF THE DIPOLES (IN WAVELENGTHS) = ');
                    dipoles_length = input(optionA_length);

           while mod(2*dipoles_length,2)~= 1
	            disp(' ');
                    disp('                       **********ERROR*********');
                    disp('   THE LENTGH OF THE DIPOLES HAS TO BE AN ODD MULTIPLE OF HALF WAVELENGTH');
                    optionA_length = sprintf('\n   LENGTH OF THE DIPOLES (IN WAVELENGTHS) = ');
                    dipoles_length = input(optionA_length);
           end;

                    optionA_horizontaldis = sprintf('   HORIZONTAL DISPLACEMENT OF THE DIPOLES ');
                    optionA_horizontaldis = sprintf('%s\n   (IN WAVELENGTHS )     ',optionA_horizontaldis);
                    optionA_horizontaldis = sprintf('%s                 = ',optionA_horizontaldis);

                   
                    horizontal_displacement = input(optionA_horizontaldis);
                    l = dipoles_length;
                    d = horizontal_displacement;

               %Mutual impedance calculations of a dipole

                   %Side-by-Side Configuration
  
                   % Based on current maximum
      
		   u0=k*d;
                   u1=k*(sqrt(d^2+l^2)+l);
                   u2=k*(sqrt(d^2+l^2)-l);

                   R21m=(ETA/(4*pi))*(2*cosint(u0)-cosint(u1)-cosint(u2));
                   X21m=-(ETA/(4*pi))*(2*sinint(u0)-sinint(u1)-sinint(u2));

                   % Based on current at the input

                   R21i = R21m/(sin(k*l/2))^2;
                   X21i = X21m/(sin(k*l/2))^2;

              if output_device == 2
                  diary(filename);
              end
                   
                   disp(' '),
                   disp('   **************PROGRAM OUTPUT***************'),
                   disp('   MUTUAL IMPEDANCE CALCULATION OF TWO DIPOLES'),
		   disp(' '),
                   disp('   SIDE-BY-SIDE CONFIGURATION'),
                   disp(' '),
                   disp([optionA_length num2str(dipoles_length,'%12.6f')]),  
                   disp('   DISTANCE BETWEEN THE TWO DIPOLES'), 
                   disp(['   (IN WAVELENGTHS)...................... = ' num2str(horizontal_displacement,'%12.6f')]), 
                   disp(' '),
                   disp(['   R21i  =     '  num2str(R21i,'%12.6f')]),
                   disp(['   X21i  =     '  num2str(X21i,'%12.6f')]),
                   disp(['   R21m  =     '  num2str(R21m,'%12.6f')]),
                   disp(['   X21m  =     '  num2str(X21m,'%12.6f')]),  

	   diary off;
 
               case {'''B''', '''b'''}
                    
                    % If user inserts option B

                    optionB_length = sprintf('\n   LENGTH OF THE DIPOLES (IN WAVELENGTHS) = ');
                    dipoles_length = input(optionB_length);

           while mod(2*dipoles_length,2)~= 1
	            disp(' ');
                    disp('                       **********ERROR*********');
                    disp('   THE LENTGH OF THE DIPOLES HAS TO BE AN ODD MULTIPLE OF HALF WAVELENGTH');
                    optionB_length = sprintf('\n   LENGTH OF THE DIPOLES (IN WAVELENGTHS) = ');
                    dipoles_length = input(optionB_length);
           end;

                    optionB_verticaldis = sprintf('   VERTICAL DISPLACEMENT OF THE DIPOLES ');
                    optionB_verticaldis = sprintf('%s\n   (IN WAVELENGTHS )     ',optionB_verticaldis);
                    optionB_verticaldis = sprintf('%s                 = ',optionB_verticaldis);

                    
                    vertical_displacement = input(optionB_verticaldis);
                    l = dipoles_length;
                    s = vertical_displacement;

                %Mutual impedance calculations of a dipole

                    %Collinear Configuration
     
                    % Based on current maximum
    
                    h=l+s;
                    v0=k*h;
                    v1=2*k*(h+l);
                    v2=2*k*(h-l);
                    v3=(h^2-l^2)/h^2;

                    R21m=-(ETA/(8*pi))*cos(v0)*(-2*cosint(2*v0)+cosint(v2)+cosint(v1)...
                         -log(v3))+(ETA/(8*pi))*sin(v0)*(2*sinint(2*v0)-sinint(v2)-sinint(v1)); 
     
                    X21m=-(ETA/(8*pi))*cos(v0)*(2*sinint(2*v0)-sinint(v2)-sinint(v1))...
                         +(ETA/(8*pi))*sin(v0)*(2*cosint(2*v0)-cosint(v2)-cosint(v1)-log(v3));

                    % Based on current at the input

                    R21i = R21m/(sin(k*l/2))^2;
                    X21i = X21m/(sin(k*l/2))^2;

         if output_device == 2
	 
                   diary(filename);
        
         end
                   disp(' '),
                   disp('   **************PROGRAM OUTPUT***************'),
                   disp('   MUTUAL IMPEDANCE CALCULATION OF TWO DIPOLES'),
		   disp(' '),
                   disp('   COLLINEAR CONFIGURATION'),
                   disp(' '),
                   disp(['   LENGTH OF THE TWO DIPOLES (IN WAVELENGTHS) = ' num2str(dipoles_length,'%12.6f')]),  
                   disp('   SEPERATION BETWEEN THE TWO DIPOLES'), 
		   disp(['   (IN WAVELENGTHS).......................... = ' num2str(vertical_displacement,'%12.6f')]), 
                   disp(' '),
                   disp(['   R21i  =     '  num2str(R21i,'%12.6f')]),
                   disp(['   X21i  =     '  num2str(X21i,'%12.6f')]),
                   disp(['   R21m  =     '  num2str(R21m,'%12.6f')]),
                   disp(['   X21m  =     '  num2str(X21m,'%12.6f')]),  
 
	   diary off;

                case {'''C''', '''c'''}
                    
                    % If user inserts option C

                    optionC_length = sprintf('\n   LENGTH OF THE TWO DIPOLES  ');
                    optionC_length = sprintf('%s\n   (IN WAVELENGTHS)              ',optionC_length);
                    optionC_length = sprintf('%s         = ',optionC_length);
                    dipoles_length = input(optionC_length);

          while mod(2*dipoles_length,2)~= 1
	            disp(' ');
                    disp('                       **********ERROR*********');
                    disp('   THE LENTGH OF THE DIPOLES HAS TO BE AN ODD MULTIPLE OF HALF WAVELENGTH');
                    optionC_length = sprintf('\n   LENGTH OF THE DIPOLES (IN WAVELENGTHS) = ');
                    dipoles_length = input(optionC_length);
           end; 

                    
                    optionC_verticaldis = sprintf('   VERTICAL DISPLACEMENT OF THE DIPOLES ');
                    optionC_verticaldis = sprintf('%s\n   (IN WAVELENGTHS )     ',optionC_verticaldis);
                    optionC_verticaldis = sprintf('%s                 = ',optionC_verticaldis);
                    optionC_horizontaldis = sprintf('   HORIZONTAL DISPLACEMENT OF THE DIPOLES ');
                    optionC_horizontaldis = sprintf('%s\n   (IN WAVELENGTHS )     ',optionC_horizontaldis);
                    optionC_horizontaldis = sprintf('%s                 = ',optionC_horizontaldis);

                    
                    vertical_displacement = input(optionC_verticaldis);
                    horizontal_displacement = input(optionC_horizontaldis);
                    l = dipoles_length;
                    d = horizontal_displacement;
                    h = vertical_displacement;

                 %Mutual impedance calculations of a dipole
 
                    % Parallel-in-Enchelon Configuration
    
                    % Based on current maximum
     
		    w0=k*h;
                    w1=k*(sqrt(d^2+h^2)+h);
                    w1p=k*(sqrt(d^2+h^2)-h);
                    w2=k*(sqrt(d^2+(h-l)^2)+(h-l));
                    w2p=k*(sqrt(d^2+(h-l)^2)-(h-l));
                    w3=k*(sqrt(d^2+(h+l)^2)+(h+l));
                    w3p=k*(sqrt(d^2+(h+l)^2)-(h+l));

                    R21m=-(ETA/(8*pi))*cos(w0)*(-2*cosint(w1)-2*cosint(w1p)+cosint(w2)...
                         +cosint(w2p)+cosint(w3)+cosint(w3p))...
                         +(ETA/(8*pi))*sin(w0)*(2*sinint(w1)-2*sinint(w1p)-sinint(w2)...
                         +sinint(w2p)-sinint(w3)+sinint(w3p));
    
                    X21m=-(ETA/(8*pi))*cos(w0)*(2*sinint(w1)+2*sinint(w1p)-sinint(w2)...
                         -sinint(w2p)-sinint(w3)-sinint(w3p))...
                         +(ETA/(8*pi))*sin(w0)*(2*cosint(w1)-2*cosint(w1p)-cosint(w2)...
                         +cosint(w2p)-cosint(w3)+cosint(w3p));

                    % Based on current at the input

                    R21i = R21m/(sin(k*l/2))^2;
                    X21i = X21m/(sin(k*l/2))^2; 

        if output_device == 2
    
                   diary(filename);
        
         end
                   disp(' '),
                   disp('   **************PROGRAM OUTPUT***************'),
                   disp('   MUTUAL IMPEDANCE CALCULATION OF TWO DIPOLES'),
		   disp(' '),
                   disp('   PARALLEL-IN-ENCHELON CONFIGURATION'),
                   disp(' '),
                   disp(['   LENGTH OF THE TWO DIPOLES (IN WAVELENGTHS) = ' num2str(dipoles_length,'%12.6f')]),  
                   disp('   DISTANCE BETWEEN THE TWO DIPOLES'),
                   disp(['   (IN WAVELENGTHS).......................... = '...
                        num2str(horizontal_displacement,'%12.6f')]),
                   disp('   RELATIVE HIGHT BETWEEN THE TWO DIPOLES'),
                   disp(['   (IN WAVELENGTHS).......................... = '...
                        num2str(vertical_displacement,'%12.6f')]), 
                   disp(' '),
                   disp(['   R21i  =     '  num2str(R21i,'%12.6f')]),
                   disp(['   X21i  =     '  num2str(X21i,'%12.6f')]),
                   disp(['   R21m  =     '  num2str(R21m,'%12.6f')]),
                   disp(['   X21m  =     '  num2str(X21m,'%12.6f')]),  
 
	   diary off;

                end    
            end
            
     
       


  

     
   


















