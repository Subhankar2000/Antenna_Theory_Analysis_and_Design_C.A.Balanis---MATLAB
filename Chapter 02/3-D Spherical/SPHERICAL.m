% ************************************************************************************
% SPHERICAL
%*************************************************************************
% This interactive MATLAB program is used to plot three-dimensional (3-D) pattern
% in spherical form.
%
% The user MUST specify the function, which represents the LINEAR FIELD pattern 
% of the function whose pattern is to be plotted, in the Matlab function
%
% f.m
%
% Also the user will be asked to specify the range of angles over which the 
% which the pattern is to be plotted:
%
% A.  0 degrees < theta < 180 degrees
% B.  0 degrees < phi < 360 degrees
%
% P. S.  Neither of the two angles need to be specified over their maximum range.
%
% The user will also have the option to remove, for visual purposes, any segment 
% of the pattern, as is demonstrated throughout the textbook. 
%
% The program can plot any of the following:
%
% 1.  Field pattern (linear scale)
% 2.  Power pattern (linear scale)
% 3.  Power pattern (in dB)
%
% Program written by:
% Bo Yang, AHE Lab 
% Department of Electrical Engineering
% Arizona State University

function Patt3d

clear all;
close all;

% Data Entry
fprintf('\n--------------------------------------------------------------------------------');
fprintf('\nThis program is a Matlab funtion that plots 3D Pattern in Spherical Coordinate');
fprintf('\n--------------------------------------------------------------------------------\n');
fprintf('\n   *** NOTICE: Please edit m-file "f.m" to input field pattern expressions!\n\n');
fprintf('   Make sure to use dot operations(".*","./" and ".^") in the expression,\n');
fprintf('   instead of regular operations("*","/" and "^")!\n\n\n');;

% Angle Increment
del_th=2;
del_phi=4;

% Plot Range
ERR = 1;
while(ERR ~= 0)
   fprintf('\nPlease specify the range of theta in the form "[theta1,theta2]"\n');
   fprintf('(in degrees, for example: [0,180])\n');
   fprintf('-----------------------------------------------------------------\n');
   thrg = str2num(input('Input: ','s'));
   s=size(thrg);
   if (s(1) ~= 1) | (s(2) ~= 2)
      fprintf('\n   *** ERROR: Please specify the range in the form "[theta1,theta2]")!\n\n','r');
   elseif (thrg(1) < 0) | (thrg(1) >= thrg(2)) | (thrg(2) > 180)
      fprintf('\n   *** ERROR: Theta must satisfys: 0 < theta1 < theta2 < 180 (degrees)")!\n\n');
   else
      ERR = 0;
   end
end

ERR = 1;
while(ERR ~= 0)
   fprintf('\nPlease specify the range of phi in the form "[phi1,phi2]"\n');
   fprintf('(in degrees, for example: [0,360])\n');
   fprintf('-----------------------------------------------------------------\n');
   phirg = str2num(input('Input: ','s'));
   s=size(phirg);
   if (s(1) ~= 1) | (s(2) ~= 2)
      fprintf('\n   *** ERROR: Please specify the range in the form "[phi1,phi2]")!\n\n','r');
   elseif (phirg(1) < 0) | (phirg(1) >= phirg(2)) | (phirg(2) > 360)
      fprintf('\n   *** ERROR: Phi must satisfy: 0 < phi1 < phi2 < 360 (degrees)")!\n\n');
   else
      ERR = 0;
   end
end

thrg(2)=round((thrg(2)-thrg(1))/del_th)*del_th+thrg(1);
theta=((thrg(1):del_th:thrg(2)))*pi/180;
phirg(2)=round((phirg(2)-phirg(1))/del_phi)*del_phi+phirg(1);
phi=(phirg(1):del_phi:phirg(2))*pi/180;

% Cut
ERR = 1;
while(ERR ~= 0)
    fprintf('\nNeed to remove part of the pattern? (Y/N)\n');
    fprintf('-------------------------------------------\n');
    ff = input('Input: ','s');
    if (ff == 'Y') | (ff == 'y')
        ERR1=1;
        while(ERR1 ~= 0)
           fprintf('\nRemove the pattern between phi-plane #1 and phi-plane #2\n'); 
           fprintf('\nPlease specify phi-plane #1 and phi-plane #2 in the form "[phi#1,phi#2]"\n');
           fprintf('(in degrees), for example: [0,90]\n');
           fprintf('\n------------------------------------------------------------------------\n');           
           phi_int = str2num(input('Input: ','s'));
           s = size(phi_int);
           if (s(1) ~= 1) | (s(2) ~= 2)
               fprintf('\nPlease specify phi-plane #1 and phi-plane #2 in the form "[phi#1 phi#2]"\n');
           elseif ((phi_int(1) < 0) | (phi_int(1) > phi_int(2)) | (phi_int(2) > 360))
               fprintf('\n   *** ERROR: Phi must satisfy: 0 < phi#1 < phi#2 < 360 (degrees)")!\n\n');
           else
               ERR1 = 0;
           end
        end
        cut = 1;
        ERR = 0;
    elseif (ff == 'N') | (ff == 'n')
        cut = 0;
        ERR = 0;
    end
end

% Linear or dB Scale
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

fprintf('Calulating ... \n\n\n');

% 3D data

if cut == 1
    phi_int=(round((phi_int-phirg(1))/del_phi)*del_phi+phirg(1))*pi/180;
    idx=find((phi>phi_int(1))&(phi<phi_int(2)));
    phi(idx)=NaN;
end
[THETA,PHI]=meshgrid(theta,phi);

% Pattern Calculation
r=f(THETA,PHI);
ratio=max(max(r));

% Scaling
r = Scale(r,ratio,scl,rmin);

% Transfer to rectangular coordinate
[X,Y,Z]=Sph2rec(THETA,PHI,r);

% Cross-section - Theta
thetac1=ones(size(theta))*theta(1);
[XCT1,YCT1,ZCT1]=CrossSecT(thetac1,phi,ratio,scl,rmin);
thetac2=ones(size(theta))*theta(length(theta));
[XCT2,YCT2,ZCT2]=CrossSecT(thetac2,phi,ratio,scl,rmin);

% Cross-section -Phi
if phi(1) ~= (phi(length(phi))-2*pi)
   phic1=ones(size(phi))*phi(1);
   [XCP1,YCP1,ZCP1]=CrossSecP(theta,phic1,ratio,scl,rmin);    
   phic2=ones(size(phi))*phi(length(phi));
   [XCP2,YCP2,ZCP2]=CrossSecP(theta,phic2,ratio,scl,rmin);    
end

% Additional Cross-section - Phi
if cut ==1
    % Plane #1
    phic3=ones(size(phi))*phi_int(1);
    [XCP3,YCP3,ZCP3]=CrossSecP(theta,phic3,ratio,scl,rmin);    
    % Plane #2
    phic4=ones(size(phi))*phi_int(2);
    [XCP4,YCP4,ZCP4]=CrossSecP(theta,phic4,ratio,scl,rmin);    
end

% 3D Plot
% Plot Pattern
h=mesh(X,Y,Z,'edgecolor','black');
set(h,'linewidth',1.1);
hold on;

% Plot Cross-section - Theta
cst1=mesh(XCT1,YCT1,ZCT1);
set(cst1,'linestyle','none');
cst2=mesh(XCT2,YCT2,ZCT2);
set(cst2,'linestyle','none');

% Plot Cross-section - Phi
if phi(1) ~= (phi(length(phi))-2*pi)
   csp1=mesh(XCP1,YCP1,ZCP1);
   set(csp1,'linestyle','none');
   csp2=mesh(XCP2,YCP2,ZCP2);
   set(csp2,'linestyle','none');
   %Boundary Line
   bl2=plot3([0 0],[0 0],[0 min(min(Z))],'k');
   % Plot Boundary Line
   Bdyplot(theta(1),phi(1),ratio,scl,rmin);
   Bdyplot(theta(1),phi(length(phi)),ratio,scl,rmin);
   Bdyplot(theta(length(theta)),phi(1),ratio,scl,rmin);
   Bdyplot(theta(length(theta)),phi(length(phi)),ratio,scl,rmin);
end

% Plot Additional Cross-section - Phi
if cut == 1
    csp3=mesh(XCP3,YCP3,ZCP3);
    set(csp3,'linestyle','none');
    csp4=mesh(XCP4,YCP4,ZCP4);
    set(csp4,'linestyle','none');
    Bdyplot(theta(1),phi_int(1),ratio,scl,rmin);
    Bdyplot(theta(1),phi_int(2),ratio,scl,rmin);
    Bdyplot(theta(length(theta)),phi_int(1),ratio,scl,rmin);
    Bdyplot(theta(length(theta)),phi_int(2),ratio,scl,rmin);
end

% Coordinate
cor=1;
% Plot Coordinate
if cor ==1
    corx=max(max(X))*1.1;
    cory=max(max(Y))*1.1;
    corz=max(max(Z))*1.1;
    corr=max([corx,cory,corz]);
    q1=plot3([0 corr],[0 0],[0 0],'k');
    q2=plot3([0 0],[0 corr],[0 0],'k');
    q3=plot3([0 0],[0 0],[0 corr],'k');
    set(q1,'linewidth',1.5);
    set(q2,'linewidth',1.5);
    set(q3,'linewidth',1.5);
    t1=text(corr*1.05,0,0,'x','horizontalalignment','center');
    t2=text(0,corr*1.05,0,'y','horizontalalignment','center');
    t3=text(0,0,corr*1.05,'z','horizontalalignment','center');
    set(t1,'fontsize',[12],'fontname','Times New Roman','fontangle','italic');
    set(t2,'fontsize',[12],'fontname','Times New Roman','fontangle','italic');
    set(t2,'fontsize',[12],'fontname','Times New Roman','fontangle','italic');
end

% Set Display
% colormap('white');
% set(gcf,'Color','white);
% shading interp;
axis off;
axis equal;
view(180-45,20);
title('Please use Matlab menu in Figure window: View -> Camera Toolbar to rotate the plot.');
fprintf('\n\nPlease use Matlab menu in Figure window: View -> Camera Toolbar to rotate the plot.\n\n');

% Transfer Spherical Coordinate to Rectangular Coordinate
function [xt,yt,zt]=Sph2rec(thetat,phit,rhot)
xt=rhot.*sin(thetat).*cos(phit);
yt=rhot.*sin(thetat).*sin(phit);
zt=rhot.*cos(thetat);

% Scaling
function fs=Scale(rs,ratios,scls,rmins)
if scls < 3
   fs=(rs/ratios).^scls;
else
   rs=(rs/ratios);
   fs=20*log10(rs);
   idxs=find(fs>rmins);
   fs(idxs)=fs(idxs)-rmins;
   idxs=find(fs<=rmins);
   fs(idxs)=0;
end

% Calulate Cross-section of the Cut - Theta
function [xc,yc,zc]=CrossSecT(thetac,phic,ratioc,sclc,rminc)
[THETAC,PHIC]=meshgrid(thetac,phic);
rc=f(THETAC,PHIC);
rc = Scale(rc,ratioc,sclc,rminc);
% Mesh the cross-section plane
[IC,JC]=size(rc);
for i=1:IC
    rc(i,JC)=0;
end
[xc,yc,zc]=Sph2rec(THETAC,PHIC,rc);

% Calulate Cross-section of the Cut - Phi
function [xc,yc,zc]=CrossSecP(thetac,phic,ratioc,sclc,rminc)
[THETAC,PHIC]=meshgrid(thetac,phic);
rc=f(THETAC,PHIC);
rc = Scale(rc,ratioc,sclc,rminc);
% Mesh the cross-section plane
[IC,JC]=size(rc);
for j=1:JC
    rc(IC,j)=0;
end
[xc,yc,zc]=Sph2rec(THETAC,PHIC,rc);

% Plot Boundary Line
function Bdyplot(thetab,phib,ratiob,sclb,rminb)
rb= Scale(f(thetab,phib),ratiob,sclb,rminb);
[xb,yb,zb]=Sph2rec(thetab,phib,rb);
plot3([0,xb],[0,yb],[0,zb],'k');

