%************************************************************************
%     ARRAYS
%************************************************************************
%    This is a MATLAB based program that computes the
%    radiation characteristics of:
%
%     I.   LINEAR ARRAYS(UNIFORM & BROADSIDE NONUNIFORM)
%     II.  PLANAR ARRAY (BROADSIDE UNIFORM)
%     II.  CIRCULAR ARRAY (BROADSIDE UNIFORM)
%
%     THE UNIFORM AND BROADSIDE NONUNIFORM LINEAR ARRAYS HAVE N ELEMENTS
%     PLACED EQUIDISTANTLY ALONG THE Z-AXIS.
%
%     BROADSIDE PLANAR UNIFORM ARRAY HAS M x N ELEMENTS PLACED
%     EQUIDISTANTLY ALONG THE X AND Y AXES  
%     
%     OPTION I.  LINEAR ARRAYS
%
%       OPTION A.  UNIFORM
%         ** CHOICES: ARRAY TYPE
%
%         1. BROADSIDE (MAXIMUM ALONG THETA = 0 DEGREES)
%         2. ORDINARY END-FIRE
%            (a). MAXIMUM ALONG THETA = 0 DEGREES
%            (b). MAXIMUM ALONG THETA = 180 DEGREES
%         3. HANSEN WOODYARD END-FIRE
%            (a). MAXIMUM ALONG THETA = 0 DEGREES
%            (b). MAXIMUM ALONG THETA = 180 DEGREES
%         4. SCANNING (MAXIMUM ALONG THETA = THETA_MAX)
%
%         ** INPUT PARAMETERS:
%
%            1. NUMBER OF ELEMENTS
%            2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%            3. DIRECTION OF ARRAY MAXIMUM (THETA_MAX IN DEGREES)
%
%         ** PROGRAM OUTPUT:
%
%            1. NORMALIZED ARRAY FACTOR
%            2. DIRECTIVITY (DIMENSIONLESS & IN dB) USING NUMERICAL 
%               INTEGRATION OF THE ARRAY FACTOR
%            3. HALF-POWER BEAMWIDTH (IN DEGREES) USING AN ITERATIVE
%               METHOD (FOR ALL MAXIMA IN THE PATTERN)
%
%       OPTION B.  NONUNIFORM BROADSIDE
%         ** CHOICES: ARRAY TYPE
%
%         1. BINOMIAL
%         2. DOLPH-TSCHEBYSCHEFF
%
%         ** BINOMIAL ARRAY INPUT PARAMETERS:
%
%         1. NUMBER OF ELEMENTS
%         2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%
%         ** DOLPH-TSCHEBYSCHEFF INPUT PARAMETERS:
%
%         1. NUMBER OF ELEMENTS
%         2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%         3. SIDE LOBE LEVEL (IN POSITIVE dB; i.e., 30 dB)
%
%         ** PROGRAM OUTPUT:
%
%         1. NORMALIZED EXCITATION COEFFICIENTS (An)
%         2. NORMALIZED ARRAY FACTOR
%         3. DIRECTIVITY (IN dB) USING NUMERICAL INTEGRATION OF THE 
%            ARRAY FACTOR
%         4. HALF-POWER BEAMWIDTH (IN DEGREES) USING AN ITERATIVE 
%            METHOD (FOR ALL MAXIMA THAT OCCUR IN THE PATTERN)
%
%     OPTION II.  PLANAR ARRAY
%
%         ** ARRAY INPUT PARAMETERS:
%
%         1. NUMBER OF ARRAY ELEMENTS IN X-DIRECTION
%         2. SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN
%            WAVELENGTHS)
%         3. NUMBER OF ARRAY ELEMENTS IN Y-DIRECTION
%         4. SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN
%            WAVELENGTHS)
%         5. MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES)
%         6. MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES)
%         7. THE PHI ANGLE (IN DEGREES) AT WHICH THE 2-D ANTENNA
%            PATTERN NEEDS TO BE EVALUATED (PHIEVAL IN DEG.)
%
%
%      NOTE: ONLY THE ELEVATION ANTENNA PATTERN IS EVALUATED.  THIS
%      ----- PATTERN RANGES FROM THETA=0 DEG. TO THETA=180 DEG. WHEREAS 
%            PHI REMAINS CONSTANT AT PHIEVAL.  IF THE PATTERN NEEDS TO
%            BE EVALUATED IN THE BACKSIDE REGION OF THE 2-D ARRAY,
%            THEN THE PROGRAM NEEDS TO BE RE-RUN FOR A NEW PHIEVAL
%            WHICH MUST BE EQUAL TO THE PREVIOUS VALUE PLUS 180 DEG.
%
%         ** PROGRAM OUTPUT:
%
%         1. NORMALIZED ARRAY FACTOR EVALUATED AT A GIVEN ANGLE 
%            (PHIEVAL)
%         2. UNIFORM PROGRESSIVE PHASE SHIFT IN X AND Y DIRECTIONS 
%            (IN DEGREES)
%         3. DIRECTIVITY (IN dB) USING NUMERICAL INTEGRATION OF THE 
%            ARRAY FACTOR
%         4. HALF-POWER BEAMWIDTH (IN DEGREES) FOR ALL MAXIMA THAT 
%            OCCUR IN THE ELEVATION PLANE OF THE 2-D ANTENNA PATTERN
%
%     OPTION III.  CIRCULAR ARRAY
%         ** CHOICES: EQUATION TYPE
%
%         1. EXPONENTIAL FUNCTION FORM
%         2. BESSEL FUNCTION FORM
%
%
%     ALL THE INPUT PARAMETERS ARE IN TERMS OF THE WAVELENGTH.
%     *****************************************************************
%
%     Written by: Seunghwan Yoon, Arizona State University,  08/12/2002
%
%     Revised by: Seunghwan Yoon, Arizona State University,  12/03/2002
%
%     *****************************************************************

function[]=ARRAYS;
close all;
clc;
option_a=0;
while((option_a~=1)&(option_a~=2)),
disp(strvcat('OUTPUT DEVICE OPTION FOR THE OUTPUT PARAMETERS',...
   'OPTION (1):SCREEN','OPTION (2):OUTPUT FILE'));   
option_a=input('OUTPUT DEVICE =');
end
if option_a==2, %OUTPUT FILE
   filename=input('INPUT THE DESIRED OUTPUT FILENAME <in single quotes> = ',...
      's');
   fid=fopen(filename,'wt');
end
disc=181;
MM=disc;
NN=disc;
theta_low=0;
theta_up=180;
phi_low=0;
phi_up=360;


% CHOICE OF LINEAR OR PLANAR ARRAY OR CIRCULAR ARRAY

option_b=0;
while ((option_b~=1)&(option_b~=2)&(option_b~=3)),
disp(strvcat('LINEAR OR PLANAR ARRAY','OPTION (1):LINEAR ARRAY','OPTION (2):PLANAR ARRAY','OPTION (3):CIRCULAR ARRAY'));
option_b=input('OPTION NUMBER =');
end;

if option_b==1, % LINEAR ARRAY  
% CHOICE OF UNIFORM OR NONUNIFORM    
option_c=0;
while ((option_c~=1)&(option_c~=2)),
disp(strvcat('UNIFORM OR NONUNIFORM ARRAY','OPTION (1):UNIFORM ARRAY','OPTION (2):NONUNIFORM ARRAY'));
option_c=input('OPTION NUMBER =');
end;
M=1800;
%MM=180;
%NN=180;
k=2*pi;
theta=linspace(0,pi,M+1);
theta3=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
phi3=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
[THETA,PHI]=meshgrid(theta3,phi3);
dtheta=pi/M;

if option_c==1, % UNIFORM ARRAY 
% CHOICE OF ARRAY TYPE
% (BROADSIDE,ORDINARY END-FIRE,HANSEN-WOODYARD END-FIRE,SCANNING)
option_d=0;
while ((option_d~=1)&(option_d~=2)&(option_d~=3)&(option_d~=4)),
disp(strvcat('ARRAY NAMES','OPTION (1):BROADSIDE ARRAY(MAXIMUM ALONG THETA = 90 DEGREES)','OPTION (2):ORDINARY END-FIRE ARRAY','OPTION (3):HANSEN-WOODYARD END-FIRE ARRAY','OPTION (4):SCANNING ARRAY'));
option_d=input('OPTION NUMBER =');
end;

if option_d==1, % BROADSIDE
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end   
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
beta=0;
psi=k.*d.*cos(theta)+beta;
psi3=k.*d.*cos(THETA)+beta;
AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);   % I used sinc function.
AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);  % sinc(x)=sin(pi*x)/(pi*x).  not, sin(x)/x.
   
elseif option_d==2, % ORDINARY END-FIRE
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
thmax=90;
while((thmax~=0)&(thmax~=180)),
   thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
end
if abs(thmax)<eps,beta=-k*d;
elseif abs(thmax-180)<eps,beta=k*d;
end
psi=k*d*cos(theta)+beta;
psi3=k*d*cos(THETA)+beta;
AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);

elseif option_d==3, % HANSEN-WOODYARD END-FIRE
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
thmax=90;
while((thmax~=0)&(thmax~=180)),
   thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
end
if abs(thmax)<eps,beta=-(k*d+pi/(Nelem-0.00001)); %what is the value(0.00001)??
elseif abs(thmax-180)<eps,beta=k*d+pi/Nelem; 
end
psi=k*d*cos(theta)+beta;
psi3=k*d*cos(THETA)+beta;
AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);

elseif option_d==4, % SCANNING
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
beta=-k*d*cos(thmax*pi/180);
psi=k*d*cos(theta)+beta;
psi3=k*d*cos(THETA)+beta;
AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);   
end;

elseif option_c==2, % NONUNIFORM ARRAY   
% CHOICE OF ARRAY TYPE
% (BINOMIAL,DOLPH-TSCHEBYSCHEFF)
option_e=0;
while ((option_e~=1)&(option_e~=2)),
   disp(strvcat('ARRAY NAMES','OPTION (1):BINOMIAL','OPTION (2):DOLPH-TSCHEBYSCHEFF'));
	option_e=input('OPTION NUMBER =');
end;

if option_e==1, % BINOMIAL
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
beta=0;
for i=1:M+1
[AF,Ncoef,Coef]=bin(theta,Nelem,d);
end
for i=1:MM+1
[AF3,Ncoef3,Coef3]=bin(THETA,Nelem,d);
end

   
elseif option_e==2, % DOLPH-TSCHEBYSCHEFF
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
beta=0;
RdB=input('SIDE LOBE LEVEL (IN dB) =');
for i=1:M+1;
[AF,Ncoef,Coef]=tscheby(theta,Nelem,d,RdB);
end
for i=1:MM+1;
[AF3,Ncoef3,Coef3]=tscheby(THETA,Nelem,d,RdB);
end

end
Coef=Coef(1:Ncoef);
Ncoef=Coef(1:Ncoef)/Coef(Ncoef);
end

U=(abs(AF)./max(abs(AF))).^2;
Prad=2*pi*sum(U.*sin(theta).*dtheta);
D=4*pi*U/Prad;
DdB=10.*log10(D+eps);
Do=max(D);
DodB=max(DdB);
if Nelem==1|max(AF)<2*min(AF)
hp=0;
thmax=0;
else
[hp,thmax]=hpbw(U,M);
end
No_maxima=length(thmax);

elseif option_b==2, % PLANAR ARRAY
%M=180;
%N=180;
k=2*pi;
dtheta=pi/MM;
dphi=2*pi/NN;

Mx=0;   
while (Mx<2),   
Mx=floor(input('NUMBER OF ELEMENTS IN THE X-DIRECTION ='));
end 
Ny=0;   
while (Ny<2),   
Ny=floor(input('NUMBER OF ELEMENTS IN THE Y-DIRECTION ='));
end 

dx=input('SPACING dx BETWEEN THE ELEMENTS (IN WAVELENGTHS) =')-1e-10;
dy=input('SPACING dy BETWEEN THE ELEMENTS (IN WAVELENGTHS) =')-1e-10;
thmax2=input('MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES - INTEGER #)=');
phimax2=input('MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES - INTEGER #)=');
phieval=input('THE PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES - INTEGER #)=');
dtor=pi/180;
betax=-k*dx*sin(dtor*thmax2)*cos(dtor*phimax2);
betay=-k*dy*sin(dtor*thmax2)*sin(dtor*phimax2);

theta=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
phi=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);

[THETA,PHI]=meshgrid(theta,phi);
AF3=af10(THETA,PHI,Mx,Ny,dx,dy,betax,betay);

Prad=sum(sum(abs(AF3).^2.*sin(THETA)*dtheta*dphi));
D1=4*pi*abs(AF3).^2/Prad;
D1dB=10.*log10(D1);
Do=4*pi*max(max(abs(AF3).^2))/Prad;
DodB=10.*log10(Do);

theta=linspace(0,pi,10*MM+1);
phi=phieval*dtor;
AF=af10(theta,phi,Mx,Ny,dx,dy,betax,betay);
D=4*pi*abs(AF).^2/Prad;
DdB=10.*log10(D);
U=(abs(AF)./max(abs(AF))).^2;

phi180=phieval*dtor-pi;
AF180=af10(theta,phi180,Mx,Ny,dx,dy,betax,betay);
D180=4*pi*abs(AF180).^2/Prad;
D180dB=10.*log10(D180);
U180=(abs(AF180)./max(abs(AF180))).^2;




elseif option_b==3, % CIRCULAR ARRAY
   option_f=0;
while ((option_f~=1)&(option_f~=2)),
disp(strvcat('FUNCTIONS','OPTION (1):EXPONENTIAL FUNCTION FORM','OPTION (2):BESSEL FUNCTION FORM'));
option_f=input('OPTION NUMBER =');
end;

%M=180;
%N=180;
k=2*pi;
dtheta=pi/MM;
dphi=2*pi/NN;

Nelem=0;   
while (Nelem<2),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end   
rad=input('RADIUS (IN WAVELENGTHS) =');
if option_f==2, % BESSEL
m=-1;   
while (m<0),   
   m=floor(input('Residuals[>2]='));
end 
end

theta0=input('MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES - INTEGER #)=');
phi0=input('MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES - INTEGER #)=');
phieval=input('THE PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES - INTEGER #)=');

if option_f==1, % EXPONENTIAL
dtor=pi/180;
theta=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
phi=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
[THETA,PHI]=meshgrid(theta,phi);
AF3=afc(THETA,PHI,theta0,phi0,Nelem,rad);
Prad=sum(sum(abs(AF3).^2.*sin(THETA)*dtheta*dphi));
D1=4*pi*abs(AF3).^2/Prad;
D1dB=10.*log10(D1);
Do=4*pi*max(max(abs(AF3).^2))/Prad;
DodB=10.*log10(Do);

theta=linspace(0,pi,10*MM+1);
phi=phieval*dtor;
AF=afc(theta,phi,theta0,phi0,Nelem,rad);
D=4*pi*abs(AF).^2/Prad;
DdB=10.*log10(D);
U=(abs(AF)./max(abs(AF))).^2;
AFdB=10.*log10(U);

phi180=phieval*dtor-pi;
AF180=afc(theta,phi180,theta0,phi0,Nelem,rad);
D180=4*pi*abs(AF180).^2/Prad;
D180dB=10.*log10(D180);
U180=(abs(AF180)./max(abs(AF180))).^2;


end

if option_f==2, % BESSEL
dtor=pi/180;
theta=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
phi=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
[THETA,PHI]=meshgrid(theta,phi);   
AF3=afc3(THETA,PHI,theta0,phi0,Nelem,rad,m);
Prad=sum(sum(abs(AF3).^2.*sin(THETA)*dtheta*dphi));
D1=4*pi*abs(AF3).^2/Prad;
D1dB=10.*log10(D1);
Do=4*pi*max(max(abs(AF3).^2))/Prad;
DodB=10.*log10(Do);

AF30=afc3(THETA,PHI,theta0,phi0,Nelem,rad,0);
Prad0=sum(sum(abs(AF30).^2.*sin(THETA)*dtheta*dphi));
D10=4*pi*abs(AF30).^2/Prad0;
D1dB0=10.*log10(D10);
Do0=4*pi*max(max(abs(AF30).^2))/Prad0;
DodB0=10.*log10(Do0);


theta=linspace(0,pi,10*MM+1);
phi=phieval*dtor;
AF=afc3(theta,phi,theta0,phi0,Nelem,rad,m);
D=4*pi*abs(AF).^2/Prad;
DdB=10.*log10(D);
U=(abs(AF)./max(abs(AF))).^2;
AFdB=10.*log10(U);

AF0=afc3(theta,phi,theta0,phi0,Nelem,rad,0);
D0=4*pi*abs(AF0).^2/Prad0;
DdB0=10.*log10(D0);
U0=(abs(AF0)./max(abs(AF0))).^2;
AFdB0=10.*log10(U0);

phi180=phieval*dtor-pi;
AF180=afc3(theta,phi180,theta0,phi0,Nelem,rad,m);
D180=4*pi*abs(AF180).^2/Prad;
D180dB=10.*log10(D180);
U180=(abs(AF180)./max(abs(AF180))).^2;


end
end
if option_a==2
   diary(filename);
end

scale=0;
while ((scale~=1)&(scale~=2)),
disp(strvcat('DIMENSIONLESS OR dB SCALE IN 3D DIRECTIVITY PLOT','OPTION (1):DIMENSIONLESS SCALE','OPTION (2):dB SCALE'));
scale=input('OPTION NUMBER =');
end;




% Let's go output!!!
disp(strvcat('********************************************************'));
disp(strvcat('PROGRAM OUTPUT'));
disp(strvcat('********************************************************'));
disp(strvcat('INPUT SPECIFICATION'));
disp(strvcat('--------------------------------------------------------'));

if option_b==1   
if option_c==1&option_d==1
disp(strvcat('UNIFORM BROADSIDE ARRAY'));
end
if option_c==1&option_d==2
disp(strvcat('UNIFORM ORDINARY END-FIRE ARRAY'));
end
if option_c==1&option_d==3
disp(strvcat('UNIFORM HANSEN-WOODYARD END-FIRE ARRAY'));
end
if option_c==1&option_d==4
disp(strvcat('UNIFORM SCANNING ARRAY'));
end
if option_c==2&option_e==1
disp(strvcat('NONUNIFORM BINOMIAL ARRAY'));
end
if option_c==2&option_e==2
disp(strvcat('NONUNIFORM DOLPH-TSCHEBYSCHEFF ARRAY'));
end
disp(['NUMBER OF ARRAY ELEMENTS = ',num2str(Nelem)]);
disp(['SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS) = ',num2str(d)]);
if option_c==1&option_d~=1
disp(['MAXIMUM NEEDS TO OCCUR AT = ',num2str(thmax)]);
end
if option_c==2&option_e==2
disp(['SIDE LOBE LEVEL (IN dB) = ',num2str(RdB)]);
end
disp('OUTPUT CHARACTERISTICS OF THE ARRAY');
disp('--------------------------------------------------------');


disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
disp(['NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
for i=1:No_maxima;
disp(['HPBW FOR MAXIMUM # =',num2str(i),'    ',num2str(hp(i)),' degrees     THMAX = ',num2str(thmax(i)),' degrees']);
end
if option_b==1&option_c==2
   if 2*round(Nelem/2)~=Nelem
      Coef(1)=2*Coef(1);
      end
disp('TOTAL EXCITATION COEFFICIENTS FOR THE ARRAY DESIGN'); 
disp(Coef);
disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO EDGE)');
Ncoef=Coef/Coef(round(Nelem/2));
disp(Ncoef);
disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO CENTER)');
Ncoef=Ncoef./Ncoef(1);
disp(Ncoef);
end
end

if option_b==2
disp(strvcat('UNIFORM PLANAR ARRAY'));
disp(['NUMBER OF ELEMENTS IN X-DIRECTION = ',num2str(Mx)]);
disp(['SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS) = ',num2str(dx)]);
disp(['NUMBER OF ELEMENTS IN Y-DIRECTION = ',num2str(Ny)]);
disp(['SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS) = ',num2str(dy)]);
disp(['MAXIMUM BEAM DIRECTION - THETA (IN DEGREES) = ',num2str(thmax2)]);
disp(['MAXIMUM BEAM DIRECTION - PHI (IN DEGREES) = ',num2str(phimax2)]);
disp(['THE 2D ANTENNA PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES) = ',num2str(phieval)]);
disp(strvcat('OUTPUT CHARACTERISTICS OF THE ARRAY'));
disp(strvcat('--------------------------------------------------------'));
disp(['PROGRESSIVE PHASE SHIFT IN X-DIRECTION = ',num2str(betax/dtor),' degrees']);
disp(['PROGRESSIVE PHASE SHIFT IN Y-DIRECTION = ',num2str(betay/dtor),' degrees']);
disp('DIRECTIVITY BASED ONLY ON THE FIELDS ABOVE THE XY-PLANE')
disp(['DIRECTIVITY = ',num2str(DodB-10*log10(2)),' dB']);
disp(['DIRECTIVITY = ',num2str(Do/2),' dimensionless']);
if max(AF3)<2*min(AF3)
hp=0;
thmax=0;
else
[hp,thmax]=hpbw(U,10*MM);
end
No_maxima=length(thmax);   
disp('DIRECTIVITY BASED ON THE FIELDS ABOVE AND BELOW THE XY-PLANE')
disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
disp(['EVALUATION PLANE: NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
for i=1:No_maxima;
disp(['HPBW FOR MAXIMUM #',num2str(i),'    ',num2str(hp(i)),' degrees     THMAX = ',num2str(thmax(i)),' degrees']);
end
end

if option_b==3
disp(strvcat('UNIFORM CIRCULAR ARRAY'));
disp(['NUMBER OF ELEMENTS = ',num2str(Nelem)]);
disp(['RADIUS (IN WAVELENGTHS) = ',num2str(rad)]);
disp(['MAXIMUM BEAM DIRECTION - THETA (IN DEGREES) = ',num2str(theta0)]);
disp(['MAXIMUM BEAM DIRECTION - PHI (IN DEGREES) = ',num2str(phi0)]);
disp(['THE 2D ANTENNA PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES) = ',num2str(phieval)]);
disp(strvcat('OUTPUT CHARACTERISTICS OF THE ARRAY'));
disp(strvcat('--------------------------------------------------------'));
disp('DIRECTIVITY BASED ONLY ON THE FIELDS ABOVE THE XY-PLANE')
disp(['DIRECTIVITY = ',num2str(DodB-10*log10(2)),' dB']);
disp(['DIRECTIVITY = ',num2str(Do/2),' dimensionless']);
if max(AF3)<2*min(AF3)
hp=0;
thmax=0;
else
[hp,theta0]=hpbw(U,10*MM);
end
No_maxima=length(theta0);   
disp('DIRECTIVITY BASED ON THE FIELDS ABOVE AND BELOW THE XY-PLANE')
disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
disp(['EVALUATION PLANE: NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
for i=1:No_maxima;
disp(['HPBW FOR MAXIMUM #',num2str(i),'    ',num2str(hp(i)),' degrees     THMAX = ',num2str(theta0(i)),' degrees']);
end
end










disp(' *** NOTE:');
disp(' THE NORMALIZED ARRAY FACTOR (in dB) IS STORED IN');
disp(' AN OUTPUT FILE CALLED  ............   ArrFac.dat');
disp(' ================================================');

diary off;

AFdB=10.*log10(U);
if option_b==1,
for i=1:901
   thetarec(i)=theta(i*2-1);
   AFdBrec(i)=AFdB(i*2-1);
   end
   for i=902:1801
   thetarec(i)=2*pi-theta(3601-i*2+3);
   AFdBrec(i)=AFdB(3601-i*2+3);

end
elseif option_b==2|option_b==3,
for i=1:MM+1
   thetarec(i)=theta(i*10-9);
   AFdBrec(i)=AFdB(i*10-9);
end
for i=MM+2:2*MM+1
   thetarec(i)=2*pi-thetarec(2*MM+2-i);
   AFdBrec(i)=AFdBrec(2*MM+2-i);
end

end

fidaf=fopen('ArrFac.dat','wt');
fprintf(fidaf,'%7.3f        %9.5f\n',[thetarec.*180/pi; AFdBrec]);
fclose(fidaf);














% PLOT THE GRAPHS
% ARRAY FACTOR

clf;
plot(theta*180/pi,AFdB,'m','linewidth',2); 
xlabel(['\theta',' (degrees)']),ylabel('ARRAY FACTOR(dB)')
grid on;
axis([0 180 max(min(AFdB)-1,-60) 1]);
t1=text(1,1,['HPBW = ',num2str(max(hp)),' (degrees)']);
set(t1,'units','normalized','position',[1 1.05],'horizontalalign','right');

if option_b==1&option_c==1&option_d==1
   s1=title('UNIFORM BROADSIDE','Fontsize',15);
   set(gca,'units','normalized');
   set(s1,'position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==2
   s2=title('UNIFORM ORDINARY END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s2,'position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==3
   s3=title('UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s3,'position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==4
   s4=title('UNIFORM SCANNING','Fontsize',15);
      set(gca,'units','normalized');
   set(s4,'position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==1
   s5=title('NONUNIFORM BINOMIAL','Fontsize',15);
      set(gca,'units','normalized');
   set(s5,'position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==2
   s6=title('NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
      set(gca,'units','normalized');
   set(s6,'position',[0 1],'horizontalalign','left');
end
if option_b==2
   s7=title('UNIFORM PLANAR','Fontsize',15);
      set(gca,'units','normalized');
   set(s7,'position',[0 1],'horizontalalign','left');
end
if option_b==3&option_f==1
   s8=title('UNIFORM CIRCULAR(EXPONENTIAL)','Fontsize',15);
      set(gca,'units','normalized');
   set(s8,'position',[0 1],'horizontalalign','left');
end
if option_b==3&option_f==2
   hold on;
   plot(theta*180/pi,AFdB0,'b--','linewidth',2);
   legend('principal+residuals','principal')
   s9=title('UNIFORM CIRCULAR(BESSEL)','Fontsize',15);
      set(gca,'units','normalized');
   set(s9,'position',[0 1],'horizontalalign','left');
end

figure;



diff=Do-min(D);
subplot(2,1,1)
plot(theta*180/pi,D,'r','linewidth',2);
xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dimensionless)')
grid on;
axis([0 180 floor(min(D)-0.1*diff-.1) ceil(Do+0.1*diff+.1)]);
t2=text(1,1,['D_0 = ',num2str(Do),' (dimensionless)']);
set(t2,'units','normalized','position',[1 1.05],'horizontalalign','right');

if option_b==1&option_c==1&option_d==1
   s1=title('UNIFORM BROADSIDE','Fontsize',15);
    set(gca,'units','normalized');
   set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==2
   s2=title('UNIFORM ORDINARY END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==3
   s3=title('UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s3,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==4
   s4=title('UNIFORM SCANNING','Fontsize',15);
      set(gca,'units','normalized');
   set(s4,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==1
   s5=title('NONUNIFORM BINOMIAL','Fontsize',15);
      set(gca,'units','normalized');
   set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==2
   s6=title('NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
      set(gca,'units','normalized');
   set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==2
   s7=title('UNIFORM PLANAR','Fontsize',15);
      set(gca,'units','normalized');
      set(s7,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==3&option_f==1
   s8=title('UNIFORM CIRCULAR(EXPONENTIAL)','Fontsize',15);
      set(gca,'units','normalized');
      set(s8,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==3&option_f==2
%         hold on;
%         plot(theta*180/pi,D0,'b--','linewidth',2); 
%         legend('principal+residuals','principal')
   s9=title('UNIFORM CIRCULAR(BESSEL)','Fontsize',15);
      set(gca,'units','normalized');
      set(s9,'units','normalized','position',[0 1],'horizontalalign','left');
end




diffdB=DodB-min(DdB);
subplot(2,1,2)
plot(theta*180/pi,DdB,'b','linewidth',2);
%if option_b==3&option_f==2
%         hold on;
%         plot(theta*180/pi,D0,'r--','linewidth',2); 
%         legend('principal+residuals','principal');
%         end
t3=text(1,1,['D_0 = ',num2str(DodB),' (dB)']);
set(t3,'units','normalized','position',[1 1.05],'horizontalalign','right');
xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dB)')
grid on;
axis([0 180 max(-50,10*floor(min(DdB)/10)) 10*ceil(DodB/10)]);



if option_b==1&option_c~=2;
   for i=1:Nelem;
   Ncoef(i)=1;
end
end
if option_b==1;
   figure;
   x=(1-Nelem)*d/2:d:(Nelem-1)*d/2;
   y=(1-Nelem)/2:1:(Nelem-1)/2;
   [AX,H1,H2]=plotyy(x,Ncoef(ceil(abs(y)+0.1)),x,beta.*y.*180./pi);
   set(get(AX(1),'Ylabel'),'String','AMPLITUDE','color','r');
set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
set(AX(1),'ycolor','r');
set(AX(2),'ycolor','b');
set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');

xlabel(['ARRAY LENGTH',' (\lambda)']);
%axis([],'color','r');
%axis([(1-Nelem)*d/2 (Nelem-1)*d/2 0 max(Ncoef)+0.1],'color','r');
grid on;
l2=legend('AMPLITUDE',2);
[hle,l3]=legend('PHASE',1); set(hle,'color',[1 1 1]);

if option_b==1&option_c==1&option_d==1
   s1=title('UNIFORM BROADSIDE','Fontsize',15);
   % set(gca,'units','normalized');
   set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==2
   s2=title('UNIFORM ORDINARY END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==3
   s3=title('UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
      set(gca,'units','normalized');
   set(s3,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==1&option_d==4
   s4=title('UNIFORM SCANNING','Fontsize',15);
      set(gca,'units','normalized');
   set(s4,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==1
   s5=title('NONUNIFORM BINOMIAL','Fontsize',15);
      set(gca,'units','normalized');
   set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_b==1&option_c==2&option_e==2
   s6=title('NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
      set(gca,'units','normalized');
   set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
end
end
if option_b==2
   
   for i=1:Mx*Ny;
   Ncoef(i)=1;

end

  figure;
   x=(1-Mx)*dx/2:dx:(Mx-1)*dx/2;
   y=(1-Mx)/2:1:(Mx-1)/2;
   subplot(2,1,1)
   [AX,H1,H2]=plotyy(x,Ncoef(ceil(abs(y)+0.1)),x,betax.*y.*180./pi);
  set(get(AX(1),'Ylabel'),'String','AMPLITUDE(X)','color','r');
  set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
  set(AX(1),'ycolor','r');
set(AX(2),'ycolor','b');
set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o'); 
xlabel(['ARRAY LENGTH ',' (\lambda)']);
%axis([(1-Mx)*dx/2 (Mx-1)*dx/2 0 max(Ncoef)+0.1]);
grid on;
l2=legend('AMPLITUDE',2);
[hle,l3]=legend('PHASE',1); set(hle,'color',[1 1 1]);

s7=title('UNIFORM PLANAR','Fontsize',15);
      set(gca,'units','normalized');
   set(s7,'units','normalized','position',[0 1],'horizontalalign','left');

x1=(1-Ny)*dy/2:dy:(Ny-1)*dy/2;
   y1=(1-Ny)/2:1:(Ny-1)/2;
   subplot(2,1,2)
   [AX,H1,H2]=plotyy(x1,Ncoef(ceil(abs(y1)+0.1)),x1,betay.*y1.*180./pi);
 set(get(AX(1),'Ylabel'),'String','AMPLITUDE(Y)','color','r');
 set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
 set(AX(1),'ycolor','r');
set(AX(2),'ycolor','b');
set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');  
xlabel(['ARRAY LENGTH ',' (\lambda)']);
%axis([(1-Ny)*dy/2 (Ny-1)*dy/2 0 max(Ncoef)+0.1]);
grid on;
l2=legend('AMPLITUDE',2);
[hle,l3]=legend('PHASE',1); set(hle,'color',[1 1 1]);

end
figure;
if option_b==1,
polar_dB(theta*180/pi,DdB,max(-60,6*floor(min(DdB)/6)),10*ceil(DodB/10),12,'-')
hold on
polar_dB(-theta*180/pi,DdB,max(-60,6*floor(min(DdB)/6)),10*ceil(DodB/10),12,'-')
title('Polar plot of Directivity (0< \phi <360 degrees)','Fontsize',15)
else
   polar_dB(theta*180/pi,DdB,max(-60,6*floor(min(DdB)/6)),10*ceil(max(DdB)/10),12,'-')
   hold on
polar_dB(-theta*180/pi,D180dB,max(-60,6*floor(min(DdB)/6)),10*ceil(max(DdB)/10),12,'-')
title(['Polar plot of Directivity (\phi= ',num2str(phieval),' degrees)'],'Fontsize',15)
   end
   
figure;
if option_b==1,
polar_dB(theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
hold on
polar_dB(-theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
title('Polar plot of Relative Directivity (0< \phi <360 degrees)','Fontsize',15)
else
   polar_dB(theta*180/pi,DdB-max(max(DdB),max(D180dB)),max(-60,6*floor(min(DdB)/6)),0,12,'-')
   hold on
polar_dB(-theta*180/pi,D180dB-max(max(DdB),max(D180dB)),max(-60,6*floor(min(DdB)/6)),0,12,'-')
title(['Polar plot of Relative Directivity (\phi= ',num2str(phieval),' degrees)'],'Fontsize',15)
   end

%Spherical Plot3D
D3=4*pi*abs(AF3).^2/Prad;

D3dB=10.*log10(D3);
D3dB=D3dB-min(min(D3dB));

%figure;
%[x,y,z]=sph2cart(PHI,pi/2-THETA,D3dB);
%surf(x,y,z);
%title('Directivity(Spherical Plot3D)','Fontsize',15);


if scale==1
   DD=D3;
   end
if scale==2
   DD=D3dB;
   end
disc=size(DD,1);
spherical_plot(DD,THETA,PHI,disc)
if scale==1
      ss=title('3D Spherical plot of Directivity (Dimensionless)','Fontsize',15);
   else
      ss=title('3D Spherical plot of Directivity (dB)','Fontsize',15);
   end
%title('3D Spherical plot of Directivity','Fontsize',15)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPBWCALC
function[hp,thmax]=hpbw(U,M)
tol=0.001;
imax=0;
j=0;
for i=1:M+1;
   if abs(U(i)-1)<tol & floor((j+2)/2)==imax+1,
      imax=imax+1;
      thmax(imax)=(i-1)/10;
   end
   if i>1 & abs(U(i)-1)<tol & U(i)>U(i-1) & j~=0,
      thmax(imax)=(i-1)/10;
   end
   if i>1,
      y(1)=U(i)-0.5;
      y(2)=U(i-1)-0.5;
      x(1)=(i-1)/10;
      x(2)=(i-2)/10;
      sign=y(1)*y(2);
      if sign<0,
         j=j+1;
         root(j)=x(2)-y(2)*(x(2)-x(1))/(y(2)-y(1));
         if j>=2 & y(2)>y(1),
            hp(imax)=root(j)-root(j-1);
         elseif j==1 & y(2)>y(1),
            hp(imax)=2.*root(j);
         end
      end
   end
end
if thmax(imax)>root(j),
   hp(imax)=2.*(180-root(j));
end

% BINOMIAL
function[AF,Ncoef,Coef]=bin(theta,Nelem,d)
if 2*floor(Nelem/2)==Nelem,Ncoef=Nelem/2;
else Ncoef=(Nelem+1)/2;
end
for i=1:Ncoef;
   Coef(i)=1;
   for j=1:Ncoef-i;
      Coef(i)=Coef(i).*(Nelem-j)./j;
   end
end
if 2*floor(Nelem/2)~=Nelem,Coef(1)=Coef(1)/2;
end
u=pi*d*cos(theta);
if 2*floor(Nelem/2)==Nelem,
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i).*cos((2.*i-1).*u);
   end
else AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i).*cos(2.*(i-1).*u);
   end
end

% TSCEBY(THETA,NELEM,D,RDB)
function[AF,Ncoef,Coef]=tscheby(theta,Nelem,d,RdB);
Ro=10^(RdB/20);
P=Nelem-1;
Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
if 2*floor(Nelem/2)==Nelem,
   Ncoef=Nelem/2;
   M=Ncoef;
for i=1:M;
   Coef(i)=0;
   for j=i:M;
      Coef(i)=Coef(i)+(-1)^(M-j)*Zo^(2*j-1)*fact(j+M-2)*(2*M-1)/(fact(j-i)*fact(j+i-1)*fact(M-j));
   end
end
elseif 2*floor((Nelem+1)/2)==Nelem+1,
   Ncoef=(Nelem+1)/2;
   M=Ncoef-1;
for i=1:M+1;
   Coef(i)=0;
   for j=i:M+1;
      if i==1,EN=2;
      else EN=1;
      end
      Coef(i)=Coef(i)+(-1)^(M-j+1)*Zo^(2*(j-1))*fact(j+M-2)*2*M/(EN*fact(j-i)*fact(j+i-2)*fact(M-j+1));
   end
end
end
u=pi*d*cos(theta);
if 2*floor(Nelem/2)==Nelem,
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i)*cos((2*i-1)*u);
   end
elseif 2*floor((Nelem+1)/2)==Nelem+1, 
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i)*cos(2*(i-1)*u);
   end
end

% FACT(IARG)   
function[f7]=fact(iarg)
f7=1;
for j=1:iarg;
f7=j*f7;
end

% PLANAR(THETA,PHI,MX,NY,DX,DY,BETAX,BETAY)
function[f10]=af10(theta,phi,Mx,Ny,dx,dy,betax,betay)
k=2*pi;
psix=k.*dx.*sin(theta).*cos(phi)+betax;
psiy=k.*dy.*sin(theta).*sin(phi)+betay;
AFx=sinc((Mx.*psix./2)./pi)./sinc((psix./2)./pi);
AFy=sinc((Ny.*psiy./2)./pi)./sinc((psiy./2)./pi);
f10=AFx.*AFy;


% CIRCULAR(THETA,PHI,THETA0,PHI0,NELEM,RAD)
function[fc]=afc(theta,phi,theta0,phi0,Nelem,rad)
k=2*pi;
dtor=pi/180;
%n=linspace(1,Nelem,Nelem);
%AF=sum(exp(i.*(k.*rad.*sin(theta).*cos(phi-phin(n))+alpha(n))));
AF=0;
for n=1:Nelem
   phin(n)=2*pi*n/Nelem;
alpha(n)=-k.*rad.*sin(dtor.*theta0).*cos(dtor.*phi0-phin(n));
AF=AF+exp(i.*(k.*rad.*sin(theta).*cos(phi-phin(n))+alpha(n)));
end
AFabs=abs(AF);
fc=AFabs;


% bessel with residuals
function[fc3]=afc3(theta,phi,theta0,phi0,Nelem,rad,m)
if theta0==0|phi0==0,
   theta0=theta0+0.000001;
   phi0=phi0-0.000001;
   end
k=2*pi;
dtor=pi/180;
rho0=rad.*sqrt((sin(theta).*cos(phi)-sin(theta0).*cos(phi0)).^2+(sin(theta).*sin(phi)-sin(theta0).*sin(phi0)).^2);
zeta=atan((sin(theta).*sin(phi)-sin(theta0).*sin(phi0))./(sin(theta).*cos(phi)-sin(theta0).*cos(phi0)));
AFpri=besselj(0,k.*rho0);
AFres=0;
for n=1:m
AFres=AFres+(i^(n.*Nelem))*besselj(n.*Nelem,k.*rho0).*cos(n.*Nelem.*zeta);
end
AFres=2*AFres;
AF=AFpri+AFres;
AFabs=abs(AF);
fc3=AFabs;


% spherical plot
function spherical_plot(r,THETA,PHI,disc)
%theta = linspace(theta_low,theta_up,disc);
%phi   = linspace(phi_low,phi_up,disc);

%[THETA,PHI] = meshgrid(theta,phi);

% spherical to rectangular conversion
x = abs(r).*sin(THETA).*cos(PHI);
y = abs(r).*sin(THETA).*sin(PHI);
z = abs(r).*cos(THETA);

% do the plot
figure; surf(x,y,z); view(135,20);
C = [.8 .8 .8]; colormap(C); axis off equal;

% Draw x, y, and z axes
set(line([1e-8;max(max(x))+3],[1e-8;1e-8],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;max(max(y))+3],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;1e-8],[1e-8;max(max(z))+3]),'Color','r');

% Label x, y, and z axes
text(max(max(x))+4,0,0,'x','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,max(max(y))+4,0,'y','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,0,max(max(z))+4,'z','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');

% Fill surface using patches
patch_1 = zeros(3,disc+1);  patch_2 = zeros(3,disc+1);
patch_1(1,1:disc) = x(1,:); patch_2(1,1:disc) = x(disc,:);
patch_1(2,1:disc) = y(1,:); patch_2(2,1:disc) = y(disc,:);
patch_1(3,1:disc) = z(1,:); patch_2(3,1:disc) = z(disc,:);
patch(patch_1(1,:),patch_1(2,:),patch_1(3,:),C);
patch(patch_2(1,:),patch_2(2,:),patch_2(3,:),C);



%-----------------------------------------------------------------------------
%       polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%
%	POLAR_DB is a MATLAB function that plots 2-D patterns in
%	polar coordinates where:
%		   0      <= THETA (in degrees) <= 360
%		-infinity <  RHO   (in dB)      <  +infinity
%
%	Input Parameters Description
%	----------------------------
%	- theta (in degrees) must be a row vector from 0 to 360 degrees
%	- rho (in dB) must be a row vector
%	- rmin (in dB) sets the minimum limit of the plot (e.g., -60 dB)
%	- rmax (in dB) sets the maximum limit of the plot (e.g.,   0 dB)
%	- rticks is the # of radial ticks (or circles) desired. (e.g., 4)
%	- linestyle is solid (e.g., '-') or dashed (e.g., '--')
%
%	Credits:
%		S. Bellofiore
%		S. Georgakopoulos
%		A. C. Polycarpou
%		C. Wangsvick
%		C. Bishop
%
%	Tabulate your data accordingly, and call polar_dB to provide the
%	2-D polar plot
%
%	Note:  This function is different from the polar.m (provided by
%	       MATLAB) because RHO is given in dB, and it can be negative
%-----------------------------------------------------------------------------

function hpol = polar_dB(theta,rho,rmin,rmax,rticks,line_style) 

% Convert degrees into radians
theta = theta * pi/180;

% Font size, font style and line width parameters
font_size  = 16;
font_name  = 'Times';
line_width = 1.5;

if nargin < 5
	error('Requires 5 or 6 input arguments.')
elseif nargin == 5 
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 1
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   font_size, ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
  % v returns the axis limits
  % changed the following line to let the y limits become negative
	hhh=plot([0 max(theta(:))],[min(rho(:)) max(rho(:))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);

% check radial limits (rticks)

 	if rticks > 5   % see if we can reduce the number
 		if rem(rticks,2) == 0
 			rticks = rticks/2;
 		elseif rem(rticks,3) == 0
 			rticks = rticks/3;
 		end
 	end

% define a circle
	th = 0:pi/50:2*pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:(length(th)-1)/4:length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);

	rinc = (rmax-rmin)/rticks;

% label r
  % change the following line so that the unit circle is not multiplied
  % by a negative number.  Ditto for the text locations.
	for i=(rmin+rinc):rinc:rmax
                is = i - rmin;
		plot(xunit*is,yunit*is,'-','color',tc,'linewidth',0.5);
		text(0,is+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
	end
% plot spokes
	th = (1:6)*2*pi/12;
	cst = cos(th); snt = sin(th);
	cs = [-cst; cst];
	sn = [-snt; snt];
	plot((rmax-rmin)*cs,(rmax-rmin)*sn,'-','color',tc,'linewidth',0.5);

% plot the ticks
	george=(rmax-rmin)/30; % Length of the ticks
        th2 = (0:36)*2*pi/72;
        cst2 = cos(th2); snt2 = sin(th2);
	cs2 = [(rmax-rmin-george)*cst2; (rmax-rmin)*cst2];
	sn2 = [(rmax-rmin-george)*snt2; (rmax-rmin)*snt2];
	plot(cs2,sn2,'-','color',tc,'linewidth',0.15); % 0.5
        plot(-cs2,-sn2,'-','color',tc,'linewidth',0.15); % 0.5


% annotate spokes in degrees
  % Changed the next line to make the spokes long enough
	rt = 1.1*(rmax-rmin);
	for i = 1:max(size(th))
		text(rt*cst(i),rt*snt(i),int2str(abs(i*30-90)),'horizontalalignment','center' );
		if i == max(size(th))
			loc = int2str(90);
		elseif i*30+90<=180
			loc = int2str(i*30+90);
                else
                        loc = int2str(180-(i*30+90-180));  
		end
		text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center' );
	end
% set viewto 2-D
	view(0,90);

% set axis limits
  % Changed the next line to scale things properly
	axis((rmax-rmin)*[-1 1 -1.1 1.1]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
  % changed the next line so negative rho are not plotted on the other side
  
  for i = 1:length(rho)
    if (rho(i) > rmin)
      if theta(i)*180/pi >=0 & theta(i)*180/pi <=90
          xx(i) = (rho(i)-rmin)*cos(pi/2-theta(i));
          yy(i) = (rho(i)-rmin)*sin(pi/2-theta(i));
      elseif theta(i)*180/pi >=90
          xx(i) = (rho(i)-rmin)*cos(-theta(i)+pi/2);
          yy(i) = (rho(i)-rmin)*sin(-theta(i)+pi/2);
      elseif theta(i)*180/pi < 0 
          xx(i) = (rho(i)-rmin)*cos(abs(theta(i))+pi/2);
          yy(i) = (rho(i)-rmin)*sin(abs(theta(i))+pi/2);
      end
    else
      xx(i) = 0;
      yy(i) = 0;
    end
  end

% plot data on top of grid
if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style);
end
if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end

