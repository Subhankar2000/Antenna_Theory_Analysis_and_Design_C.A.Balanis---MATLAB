%*******************************************************************************
%      PYRAMIDAL HORN: ANALYSIS
%**************************************************************************
%
%     THIS PROGRAM COMPUTES FOR A PYRAMIDAL HORN THE:
%
%     I.    FAR-FIELD E- AND H-PLANE AMPLITUDE PATTERNS* BASED ON THE
%           THEORY OF SECTION 13.4 EQUATIONS (13-46) - (13-48c)
%     II.   DIRECTIVITY (IN dB) BASED ON EQUATION (13-52)
%     III.  DIRECTIVITY (IN dB) OF THE CORRESPONDING E-PLANE SECTORAL
%           HORN BASED ON EQUATION (13-19)
%     IV.   DIRECTIVITY (IN dB) OF THE CORRESPONDING H-PLANE SECTORAL
%           HORN BASED ON EQUATION (13-41)
%
%        ** INPUT PARAMETERS:
%
%       1.  RHO1 (IN WAVELENGTHS)
%       2.  RHO2 (IN WAVELENGTHS)
%       3.  WAVEGUIDE DIMENSIONS a & b (IN WAVELENGTHS)
%       4.  HORN APERTURE DIMENSIONS a1 & b1 (IN WAVELENGTHS)
%
%
%        ** NOTE: REFER TO FIGURE 13.18 FOR THE GEOMETRY.
%                 THE E- AND H-PLANE AMPLITUDE PATTERNS ARE STORED IN
%                 TWO DATA FILES NAMELY E-Plane Horn.dat AND H-Plane Horn.dat,
%                 RESPECTIVELY.
%**************************************************************************
%     Written by: Keith B. Washington, Arizona State University
%
%**************************************************************************

function []=horn;

close all;
% Input parameters

disp(' ')
disp('E-Plane and H-Plane Horn Specifications');
disp('--------------------------------------------')

R1=[]; R2=[];
while isempty(R1)|~isa(R1,'double')|R1<=0,
   R1 = input('rho1(in wavelengths) = ');
end;

while isempty(R2)|~isa(R2,'double')|R2<=0,
   R2 = input('rho2(in wavelengths) = ');
end;

disp(' ');
disp('Waveguide Dimensions');
disp('--------------------------------------------');

a=[]; b=[];
while isempty(a)|~isa(a,'double')|a<=0,
   a = input('a(in wavelengths) = ');
end;

while isempty(b)|~isa(b,'double')|b<=0,
   b = input('b(in wavelengths) = ');
end;  

disp(' ');
disp('Horn Aperture Dimensions')
disp('--------------------------------------------')

a1=[]; b1=[];
while isempty(a1)|~isa(a1,'double')|a1<=0,
   a1 = input('a1(in wavelengths) = ');
end;

while isempty(b1)|~isa(b1,'double')|b1<=0,
   b1 = input('b1(in wavelengths) = ');
end;

disp(' ')

opt=[];
while isempty(opt)|(opt~=1&opt~=2),
   opt=input(['Select output option\n','--------------------------------------------\n', ...
              '1. Screen\n','2. Output file\n','->']);
end;

if opt==2,
   fn = input('Input the Desired Filename:\n','s');
end;

disp(' ')
disp('please wait')


% compute internal functions

u = (1/sqrt(2))*((sqrt(R2)/a1)+(a1/sqrt(R2)));

v = (1/sqrt(2))*((sqrt(R2)/a1)-(a1/sqrt(R2)));

u = Fresnel(u);

v = Fresnel(v);

DH = 4*pi*b*R2/a1*((real(u)-real(v))^2 + (imag(u)-imag(v))^2);

w = Fresnel(b1/sqrt(2*R1));

DE = 64*a*R1/(pi*b1)*((real(w))^2 + (imag(w))^2);

DP = pi/(32*a*b)*DE*DH;

k = 2*pi;

Emax = 0;

Hmax = 0;


% E and H plane Outputs


% E-Plane Amplitude

fid1 = fopen('E-Plane Horn.dat','wt');

fprintf(fid1,'                         E-Plane Horn\n');
fprintf(fid1,'Theta     E-Theta\n');


for(theta = 0:0.5:360);

    I = theta*2 + 1;

    theta = theta*pi/180;

    phi = pi/2;

    ky = k*sin(theta);

    kxp = pi/a1;

    kxdp = -pi/a1;

    t1 = sqrt(1/(pi*k*R1))*(-k*b1/2-ky*R1);

    t2 = sqrt(1/(pi*k*R1))*(k*b1/2-ky*R1);

    t1p = sqrt(1/(pi*k*R2))*(-k*a1/2-pi/a1*R2);

    t2p = sqrt(1/(pi*k*R2))*(k*a1/2-pi/a1*R2);

    t1dp = -t2p;

    t2dp = -t1p;

    I1 = .5*sqrt(pi*R2/k)*(exp(j*R2/(2*k)*kxp^2)*(Fresnel(t2p) - Fresnel(t1p)) + exp(j*R2/(2*k)*kxdp^2)*(Fresnel(t2dp) - Fresnel(t1dp)));

    I2 = sqrt(pi*R1/k) * exp(j*R1/(2*k)*ky^2) * (Fresnel(t2) - Fresnel(t1));

    y(I) = (1 + cos(theta))*I1*I2;

    y(I) = abs(y(I));

end

for(I = 1:721)
    if(y(I) > Emax)
       Emax = y(I);
    end
end

for(I = 1:721)
    if(y(I) <= 0)
       Edb = -100;
    else
       Edb = 20*log10(abs(y(I))/Emax);
    end

     theta = (I-1)/2;
     x(I)=theta;
     q1(I)=Edb;
    fprintf(fid1,'%6.1f %12.6f\n',theta,Edb);

end

fclose(fid1);

%subplot(2,1,1)
%plot(x,q);
%xlabel('Theta (degrees)');
%ylabel('Field Pattern (dB)');
%title('E-Plane');



% H-Plane Amplitude

fid2 = fopen('H-Plane Horn.dat','wt');

fprintf(fid2,'                         H-Plane Horn\n');
fprintf(fid2,'Theta     E-Phi\n');

for(theta = 0:0.5:360);

    I = theta*2 + 1;

    theta = theta*pi/180;

    phi = 0;

    kxp = k*sin(theta) + pi/a1;

    kxdp = k*sin(theta) - pi/a1;

    t1 = sqrt(1/(pi*k*R1))*(-k*b1/2);

    t2 = sqrt(1/(pi*k*R1))*(k*b1/2);

    t1p = sqrt(1/(pi*k*R2))*(-k*a1/2-kxp*R2);

    t2p = sqrt(1/(pi*k*R2))*(k*a1/2-kxp*R2);

    t1dp = sqrt(1/(pi*k*R2))*(-k*a1/2-kxdp*R2);

    t2dp = sqrt(1/(pi*k*R2))*(k*a1/2-kxdp*R2);

    I1 = .5*sqrt(pi*R2/k)*(exp(j*R2/(2*k)*kxp^2)*(Fresnel(t2p) - Fresnel(t1p)) + exp(j*R2/(2*k)*kxdp^2)*(Fresnel(t2dp) - Fresnel(t1dp)));

    I2 = sqrt(pi*R1/k) * exp(j*R1/(2*k)*ky^2) * (Fresnel(t2) - Fresnel(t1));

    y(I) = (1 + cos(theta))*I1*I2;

    y(I) = abs(y(I));

end

for(I = 1:721)
    if(y(I) > Hmax)
       Hmax = y(I);
    end
end

for(I = 1:721)
    if(y(I) <= 0)
        Hdb = -100;
    else
        Hdb = 20*log10(abs(y(I))/Hmax);
     end

     theta = (I-1)/2;
     x(I)=theta;
     q2(I)=Hdb;
    fprintf(fid2,'%6.1f %12.6f\n',theta,Hdb);

end

fclose(fid2);

% Figure 1
% ********
ha=plot(x,q1); set(ha,'linestyle','-','linewidth',2);
hold on; hb=plot(x,q2,'r--'); set(hb,'linewidth',2);
xlabel('Theta (degrees)');
ylabel('Field Pattern (dB)');
title('Horn Analysis');
legend('E-Plane','H-Plane');
grid on;
axis([0 360 -60 0]);

% Figure 2
% *********
figure(2);
ht1=elevation(x*pi/180,q1,-60,0,4,'b-');
hold on;
ht2=elevation(x*pi/180,q2,-60,0,4,'r--');
set([ht1 ht2],'linewidth',2);

legend([ht1 ht2],{'E-plane','H-plane'});
title('Field patterns');


% Directivity Output

if exist('fn'),
   fid3 = [1 fopen(fn,'wt')];
else
   fid3=1;
end;

for fin=1:length(fid3),
   fprintf(fid3(fin),'*******************************************************************************\n');
   fprintf(fid3(fin),'PROGRAM OUTPUT\n');
   fprintf(fid3(fin),'*******************************************************************************\n\n');

   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'Pyramidal Horn\n');
   fprintf(fid3(fin),'--------------\n');
   fprintf(fid3(fin),'Directivity = %.2f dB\n',10*log10(DP));
   fprintf(fid3(fin),'Directivity = %.2f dimensionless\n\n',DP);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'E-Plane Sectoral Horn\n');
   fprintf(fid3(fin),'---------------------\n');
   fprintf(fid3(fin),['Directivity = ',num2str(10*log10(DE)),' dB\n']);
   fprintf(fid3(fin),['Directivity = ',num2str(DE),' dimensionless\n']);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'H-Plane Sectoral Horn\n');
   fprintf(fid3(fin),'---------------------\n');
   fprintf(fid3(fin),['Directivity = ',num2str(10*log10(DH)),' dB\n']);
   fprintf(fid3(fin),['Directivity = ',num2str(DH),' dimensionless\n']);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'*** Note:\n');
   fprintf(fid3(fin),'    The E-Plane Amplitude Pattern is stored in E-Plane Horn.dat\n');
   fprintf(fid3(fin),'    The H-Plane Amplitude Pattern is stored in H-Plane Horn.dat\n');
end;

fclose('all');


% Fresnel Subfunction

function[y] = Fresnel(x);

A(1) = 1.595769140;
A(2) = -0.000001702;
A(3) = -6.808508854;
A(4) = -0.000576361;
A(5) = 6.920691902;
A(6) = -0.016898657;
A(7) = -3.050485660;
A(8) = -0.075752419;
A(9) = 0.850663781;
A(10) = -0.025639041;
A(11) = -0.150230960;
A(12) = 0.034404779;

B(1) = -0.000000033;
B(2) = 4.255387524;
B(3) = -0.000092810;
B(4) = -7.780020400;
B(5) = -0.009520895;
B(6) = 5.075161298;
B(7) = -0.138341947;
B(8) = -1.363729124;
B(9) = -0.403349276;
B(10) = 0.702222016;
B(11) = -0.216195929;
B(12) = 0.019547031;

CC(1) = 0;
CC(2) = -0.024933975;
CC(3) = 0.000003936;
CC(4) = 0.005770956;
CC(5) = 0.000689892;
CC(6) = -0.009497136;
CC(7) = 0.011948809;
CC(8) = -0.006748873;
CC(9) = 0.000246420;
CC(10) = 0.002102967;
CC(11) = -0.001217930;
CC(12) = 0.000233939;

D(1) = 0.199471140;
D(2) = 0.000000023;
D(3) = -0.009351341;
D(4) = 0.000023006;
D(5) = 0.004851466;
D(6) = 0.001903218;
D(7) = -0.017122914;
D(8) = 0.029064067;
D(9) = -0.027928955;
D(10) = 0.016497308;
D(11) = -0.005598515;
D(12) = 0.000838386;


if(x==0)
   y=0;
   return
elseif(x<0)
   x=abs(x);
   x=(pi/2)*x^2;
   F=0;
   if(x<4)
      for(k=1:12)
         F=F+(A(k)+j*B(k))*(x/4)^(k-1);
      end
      y = F*sqrt(x/4)*exp(-j*x);
      y = -y;
      return
   else
      for(k=1:12)
         F=F+(CC(k)+j*D(k))*(4/x)^(k-1);
      end
      y = F*sqrt(4/x)*exp(-j*x)+(1-j)/2;
      y =-y;
      return
   end
else
   x=(pi/2)*x^2;
   F=0;
   if(x<4)
      for(k=1:12)
         F=F+(A(k)+j*B(k))*(x/4)^(k-1);
      end
      y = F*sqrt(x/4)*exp(-j*x);
      return
   else
      for(k=1:12)
         F=F+(CC(k)+j*D(k))*(4/x)^(k-1);
      end
      y = F*sqrt(4/x)*exp(-j*x)+(1-j)/2;
      return
   end
end
