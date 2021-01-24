function []=taylor_cheby;

% this program synthesizes a line source using the Taylor method based on the 
% Chebyshev error

% Written by Marios Gkatzianas
% Arizona State University, September 2002

close all;
clear all;

nbar=0;
while (nbar<=1),
   if (nbar==1),
      disp('nbar must be greater than 1');
   end;   
   nbar=input('Specify nbar\n');
end;

side_lobe=0;
while (side_lobe<=0),
   side_lobe=input([ 'Specify desired side lobe level for the first ',... 
                    num2str(nbar),' lobes (in dB)\n']);
end;   

R0=10^(side_lobe/20);  % voltage ratio
A=1/pi*acosh(R0);
sig=nbar/sqrt(A^2+(nbar-0.5)^2);

len=0;
while (len<=0),
   len=input('Specify the length of the line source (in wavelengths)\n');
end;   

N=0;
while (N<=0),
   N=input('Select number of sample points for line source\n');   
end;
zs=linspace(-len/2,len/2,N);

n=1:(nbar-1);
u_con=pi*sig*sqrt(A^2+(n-0.5).^2);

if (len>=nbar),
   u_decr=(nbar:floor(len))*pi;
end;   

u=[u_con u_decr];
theta_null=acos(u/(pi*len))*180/pi; 
theta_null_sort=sort([theta_null 180-theta_null]);

p=1:(nbar-1);
c1=(fact(nbar-1))^2./(fact(nbar-1+p).*fact(nbar-1-p));
temp=ones(size(p));

for m=1:(nbar-1),
   temp=temp.*(1-(p*pi/u(m)).^2);
end;   
SF_coef=c1.*temp;

Is=ones(size(zs));

Ntheta=0;
while (Ntheta<=0),
   Ntheta=input('Give number of sample points for the angle domain [0 180]\n');
end;   
theta=linspace(0,pi,Ntheta);

u2=pi*len*cos(theta);
SF=sinc(u2/pi);

for q=1:(nbar-1),
   Is=Is+2*SF_coef(q)*cos(2*pi*q*zs/len);
   SF=SF.* (1-(u2/u(q)).^2)./(1-(u2/(q*pi)).^2);
end;
Is=Is/len;

thres=0;
while (thres>=0),   
   thres=input('Give threshold for SF (in -dB)\n');
end;   

SF_db=20*log10(abs(SF)/max(abs(SF)));
SF_db(SF_db<=thres)=thres;

% Figure 1
% ********
stem(p,abs(SF_coef));
set(gca,'xtick',[1:(nbar-1)]);

% Figure 2
% ********
figure(2);
plot(zs,abs(Is)/max(abs(Is))); grid;
title('Normalized current distribution');
set(gca,'xlim',[-len/2 len/2]); %,'xtick',linspace(-len/2,len/2,6));
xlabel('z/\lambda','fontsize',12);

% Figure 3
% ********
figure(3);
plot(theta*180/pi,SF_db); grid; hold on;
xlabel('\theta (in \circ)'); ylabel('SF (dB)');
title('Normalized space factor');
h=line([0 180],[-side_lobe -side_lobe]); set(h,'linestyle','--','color','r');

disp('The SF nulls (in degrees) are:');
disp(num2str(theta_null_sort','%6.2f'));


function [y]=fact(x);

for n=1:length(x),
   y(n)=prod(1:x(n));
end;   