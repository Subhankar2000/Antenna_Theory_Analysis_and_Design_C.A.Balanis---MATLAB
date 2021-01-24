%**********************************************************************
%  Reflector.m
%********************************************************************
%  The program computes the radiation by a parabolic reflector.
%  It can be used to investigate the:
%   *   Directivity
%   *   Beamwidth
%   *   Sidelobe Level
%   *   Polarization
%   *   Near-to-Far Zone Fields
%     *****************************************************************
%
%   TYPICAL INPUT DATA                                      Values                 
%   ...................................................................
%   Frequency (GHz)                                             3
%   ...................................................................
%   REFLECTOR INFORMATION:
%   ...................................................................
%   Focus (meters)                                              6
%   Diameter (meters)                                           6
%   Offset height (meters)                                      0
%   ....................................................................
%   FEED INFORMATION:
%   ....................................................................
%   Feed location of x axis (meters)                            0
%   Feed location of y axis (meters)                            0
%   Feed location of z axis (meters)                            6
%   Feed tilt angle (degrees)                                   0
%   Number of integration points per wavelength (typically 4)   4
%   q parameter for feed pattern [f(theta)=cos(theta)^q]        2.3
%   .....................................................................
%   OBSERVATION INFORMATION
%   .....................................................................
%   Radial distance to observation point (typically 1,000 m)    1000
%   Minimum theta angle (degrees)                               0
%   Maximum theta angle (degrees)                               30
%   Phi angle (degrees)                                         0
%   Number of observation angles                                31
%   [typically =theta(max)-theta(min)+1 for one-degree increments] 
%   .....................................................................
%   POLARIZATION INFORMATION
%   .....................................................................
%   X POL                                                       1
%   Y POL                                                       1
%   Phase difference PSI between X POL and Y POL (degrees)      90
%   .....................................................................
%************************************************************************
%     Written by: Seunghwan Yoon, Arizona State University,  12/06/2002
%
%************************************************************************
function[]=reflector;
close all;
clc;
reflant=imread('reflector.bmp');
imagesc(reflant); colormap(gray); 
axis off;
disp(strvcat('********************************************************************'));
disp(strvcat('Input parameters                                  Values'));
disp(strvcat('********************************************************************'));
freq=  input(' Frequency(GHz) =                                              ');
disp(strvcat('--------------------------------------------------------------------'));
disp(strvcat('Reflector information'));
disp(strvcat('--------------------------------------------------------------------'));
foc=   input(' Focus(meters) =                                               ');
dia=   input(' Diameter(meters) =                                            ');
xoset= input(' Offset height(meters) =                                       ');
disp(strvcat('--------------------------------------------------------------------'));
disp(strvcat('Feeder information'));
disp(strvcat('--------------------------------------------------------------------'));
xf=    input(' Feed location of x axis(meters) =                             ');
yf=    input(' Feed location of y axis(meters) =                             ');
zf=    input(' Feed location of z axis(meters) =                             ');
fdta=  input(' Feed tilt angle(degrees) =                                    ');
ndel=  input(' Number of integration points per wavelength(typically 4)=     ');
qq=    input(' q parameter for feed pattern[f(theta)=cos(theta)^q] =         ');
disp(strvcat('--------------------------------------------------------------------'));
disp(strvcat('Observation information '));
disp(strvcat('--------------------------------------------------------------------'));
disp(strvcat(' Radial distance to observation'));
rob=   input(' (typically 1000 meters for far field) =                       ');
thstart=input(' Minimum theta angle(degrees) =                                ');
thetaof=input(' Maximum theta angle(degrees) =                                ');
phiob= input(' Phi angle(degrees) =                                          ');
disp(strvcat(' Number of observation angle'));
nob=   input(' (typically=theta(max)-theta(min)+1 for 1 degree increments) = ');
disp(strvcat('--------------------------------------------------------------------'));
disp(strvcat('Polarization information'));
disp(strvcat('--------------------------------------------------------------------'));
apol=  input(' X POL =                                                       ');
bpol=  input(' Y POL =                                                       ');
psi=   input(' Phase difference PSI between X POL and Y POL(degrees) =       ');
disp(strvcat('********************************************************************'));

freq=freq*1e+9;

%degree to radian
rad=pi/180;
fdta=fdta*rad;
thetaof=thetaof*rad;
thstart=thstart*rad;
phiob=phiob*rad;
psi=psi*rad;


c0=3*1e+8;                      %the speed of electromagnetic wave in free space
rlam=c0/freq;                   %the wavelength in free space
del=rlam/ndel;

xfac=-1/2/foc;
mx=floor(dia/del);
my=mx;

%input power
qe=qq;
qh=qq;


xxx=sqrt(apol*apol+bpol*bpol);
apol=apol/xxx;
bpol=bpol/xxx;
hwk=exp(j*psi);
hwkn=conj(hwk);
hwk1=exp(j*(psi+pi));
hwkn1=conj(hwk1);

pwrin=(qe+qh+1)/(60*(2*qe+1)*(2*qh+1));

xfactor=sqrt(pwrin*30);

%Loop
h = waitbar(0,'Program is running. Please wait...');
for i=1:nob,
   thetaob=(i-1)*(thetaof-thstart)/(nob-1)+thstart;
   ctob=cos(thetaob);
   stob=sin(thetaob);
   cpob=cos(phiob);
   spob=sin(phiob);
   
   x3=rob*stob*cpob;
   y3=rob*stob*spob;
   z3=rob*ctob;
   
   hsumx=0;
   hsumy=0;
   hsumz=0;
   icount=0;
   
% Integration   
   for jxx=1:mx+1,
      x2=-dia/2+(jxx-1)*del;
      for jyy=1:my+1,
         y2=-dia/2+(jyy-1)*del+xoset;
         if (x2^2+(y2-xoset)^2)<(dia/2)^2
         
        
         z2=((x2^2+y2^2)/4/foc);
         
         y2m=y2+del/2;
         z2m=((x2^2+y2m^2)/4/foc);
         
         yphn=sqrt((del/2)^2+(z2m-z2)^2);
         yph(1)=0;
         yph(2)=del/2/yphn;
         yph(3)=(z2m-z2)/yphn;
         
         nnorm=sqrt(1+(x2^2+y2^2)*(xfac^2));
         norm(1)=x2*xfac/nnorm;
         norm(2)=y2*xfac/nnorm;
         norm(3)=1/nnorm;
         
         [xph]=vector1(yph,norm);
         
         icount=icount+1;
         
         [hx,hy,hz]=incfld(x2,y2,z2,rlam,foc,xf,yf,zf,fdta,qq,apol,bpol,psi);
         
         [jx,jy,jz]=vector(hx,hy,hz,norm(1),norm(2),norm(3));
         
         hj(1)=jx;
         hj(2)=jy;
         hj(3)=jz;
         
         [hxp]=hdot(xph,hj);

         [hyp]=hdot(yph,hj);

         [hzp]=hdot(norm,hj);
         
         p23(1)=x3-x2;
         p23(2)=y3-y2;
         p23(3)=z3-z2;
         
         dist1=sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2);
         
         [pp23(1)]=dot(xph,p23);
         
         [pp23(2)]=dot(yph,p23);
         
         [pp23(3)]=dot(norm,p23);

         if pp23(3)>=dist1
            dist1=pp23(3)+1;
         end
         
         theta2=acos(pp23(3)/(dist1+(1e-8)));
         phi2=atan2(pp23(2),pp23(1)+(1e-8));
                 
         [hex,hey,hez]=patchie(dist1,theta2,phi2,hxp,hyp,hzp,del,del,rlam);
         
         hx=hex*xph(1)+hey*yph(1)+hez*norm(1);
         hy=hex*xph(2)+hey*yph(2)+hez*norm(2);
         hz=hex*xph(3)+hey*yph(3)+hez*norm(3);
         
         hsumx=hsumx+hx;
         hsumy=hsumy+hy;
         hsumz=hsumz+hz;
      else
      end
      
      end
   end
   
   
   hsumth=ctob*cpob*hsumx+ctob*spob*hsumy-stob*hsumz;
   
   hsumphi=-spob*hsumx+cpob*hsumy;
   
   hfact11=apol*hwk1*cpob+bpol*spob;
   hfact22=-apol*hwk1*spob+bpol*cpob;
   hfact33=apol*hwkn1*spob-bpol*cpob;
   hfact44=apol*hwkn1*cpob+bpol*spob;
   
   href=hsumth*conj(hfact11)+hsumphi*conj(hfact22);
   hcr=hsumth*conj(hfact33)+hsumphi*conj(hfact44);
   
   href=href*rob/xfactor/(16*pi*pi)*4*pi;
   hcr=hcr*rob/xfactor/(16*pi*pi)*4*pi;
   ex=ctob*cpob*hsumth-spob*hsumphi;
   thetaobi(i)=thetaob;
   hrefi(i)=20*log10(abs(href));
   hcri(i)=20*log10(abs(hcr));
   waitbar(i/nob,h)
   end
close(h)
 
patt = [thetaobi/rad; hrefi; hcri];
fid = fopen('Reflector_Pattern.dat','w');
fprintf(fid,'THETA(degrees)    CO-POL(dB)   CROSS-POL(dB)\n');
fprintf(fid,'%6.2f                      %6.3f             %6.3f\n',patt);
fclose(fid);   
dmaxdb=max(hrefi);
dmaxdim=10^(dmaxdb/10);
disp('************************** Program is done **************************');
disp('The Maximum Directivity Is:');
disp([num2str(dmaxdim),' (dimensionless); ',num2str(dmaxdb),' (dB)']);
disp(' ');
disp('=====================================================================');
disp('NOTE:');
disp('The Directivity Pattern(in dB) Is Stored In');
disp('An Output File Called  ....... Reflector_Pattern.dat');
disp('=====================================================================');

%Plot the graphs
figure;
%subplot(2,2,1)
plot(thetaobi/rad,hrefi,'r','linewidth',2)
title('Parabolic Reflector (CO-POL)')
xlabel(['\theta',' (degrees)']),ylabel('Directivity(dB)')
grid on;
axis([thstart/rad thetaof/rad 10*floor(min(hrefi)/10) 10*ceil(max(hrefi)/10)]);

figure;
%subplot(2,2,2)
plot(thetaobi/rad,hcri,'b','linewidth',2)
title('Parabolic Reflector (CROSS-POL)')
xlabel(['\theta',' (degrees)']),ylabel('Directivity(dB)')
grid on;
axis([thstart/rad thetaof/rad 10*floor(min(hcri)/10) 10*ceil(max(hcri)/10)]);

figure;
%subplot(2,2,3)
polar_dB(thetaobi/rad,hrefi,10*floor(min(hrefi)/10),10*ceil(max(hrefi)/10),...
   (10*ceil(max(hrefi)/10)-10*floor(min(hrefi)/10))/10,'r-')
hold on;
polar_dB(-thetaobi/rad,hrefi,10*floor(min(hrefi)/10),10*ceil(max(hrefi)/10),...
   (10*ceil(max(hrefi)/10)-10*floor(min(hrefi)/10))/10,'r-')
title('Absolute Pattern of Parabolic Reflector (CO-POL)')

% Normalized
figure;
%subplot(2,2,3)
nhrefi=hrefi-max(hrefi);
polar_dB(thetaobi/rad,nhrefi,10*floor(min(nhrefi)/10),10*ceil(max(nhrefi)/10),...
   (10*ceil(max(nhrefi)/10)-10*floor(min(nhrefi)/10))/10,'r-')
hold on;
polar_dB(-thetaobi/rad,nhrefi,10*floor(min(nhrefi)/10),10*ceil(max(nhrefi)/10),...
   (10*ceil(max(nhrefi)/10)-10*floor(min(nhrefi)/10))/10,'r-')
title('Normalized Pattern of Parabolic Reflector (CO-POL)');

figure;
%subplot(2,2,4)
polar_dB(thetaobi/rad,hcri,10*floor(min(hcri)/10),10*ceil(max(hcri)/10),...
   (10*ceil(max(hcri)/10)-10*floor(min(hcri)/10))/10,'b-');
hold on;
polar_dB(-thetaobi/rad,hcri,10*floor(min(hcri)/10),10*ceil(max(hcri)/10),...
   (10*ceil(max(hcri)/10)-10*floor(min(hcri)/10))/10,'b-');
title('Pattern of Parabolic Reflector (CROSS-POL)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%incfld
function[hx,hy,hz]=incfld(x2,y2,z2,rlam,foc,xf,yf,zf,fdta,qq,apol,bpol,psi)
qe=qq;
qh=qq;
ca=cos(fdta);
sa=sin(fdta);
x2=-x2;
xf=-xf;
xp=x2-xf;
yp=y2-yf;
zp=zf-z2;
xxf=xp;
yyf=ca*yp-sa*zp;
zzf=sa*yp+ca*zp;
r=sqrt(xxf*xxf+yyf*yyf+zzf*zzf);
thetap=acos(zzf/r);
phip=atan2(yyf,xxf+1e-8);
sp=sin(phip);
cp=cos(phip);
st=sin(thetap);
ct=cos(thetap);
hh=exp(-j*2*pi*r/rlam)/r;
hwk=apol*exp(j*psi);
ve=abs(cos(thetap))^qe;
vh=abs(cos(thetap))^qh;
heth=hh*(hwk*cp+bpol*sp)*ve;
heph=hh*(-hwk*sp+bpol*cp)*vh;

hrp=0;
htp=-heph/120/pi;
hpp=heth/120/pi;
hxp=st*cp*hrp+ct*cp*htp-sp*hpp;
hyp=st*sp*hrp+ct*sp*htp+cp*hpp;
hzp=ct*hrp-st*htp;

hhyp=ca*hyp+sa*hzp;
hhzp=-sa*hyp+ca*hzp;

hhxp=-hxp;
hhzp=-hhzp;

hx=hhxp;
hy=hhyp;
hz=hhzp;




%vector
function[jx,jy,jz]=vector(hx,hy,hz,nx,ny,nz)
jx=2*(ny*hz-nz*hy);
jy=2*(nz*hx-nx*hz);
jz=2*(nx*hy-ny*hx);

%vector1
function[hvec3]=vector1(hvec1,hvec2)
hvec3(1)=hvec1(2)*hvec2(3)-hvec1(3)*hvec2(2);
hvec3(2)=hvec1(3)*hvec2(1)-hvec1(1)*hvec2(3);
hvec3(3)=hvec1(1)*hvec2(2)-hvec1(2)*hvec2(1);

%vectorh
function[hvec3]=vectorh(hvec1,hvec2)
hvec3(1)=hvec1(2)*hvec2(3)-hvec1(3)*hvec2(2);
hvec3(2)=hvec1(3)*hvec2(1)-hvec1(1)*hvec2(3);
hvec3(3)=hvec1(1)*hvec2(2)-hvec1(2)*hvec2(1);

%hdot
function[hscalar]=hdot(hvec1,hvec2)
hscalar=hvec1(1)*hvec2(1)+hvec1(2)*hvec2(2)+hvec1(3)*hvec2(3);

%dot
function[hscalar]=dot(hvec1,hvec2)
hscalar=hvec1(1)*hvec2(1)+hvec1(2)*hvec2(2)+hvec1(3)*hvec2(3);

%patch
function[hex,hey,hez]=patchie(dist1,theta2,phi2,hxp,hyp,hzp,xd,yd,rlam)
knum=2*pi/rlam;
zz0=120*pi;
ct=cos(theta2);
st=sin(theta2);
sp=sin(phi2);
cp=cos(phi2);
hh=exp(-j*2*pi*dist1/rlam)/dist1;
hbpol=hyp;
hwk=hxp;

xx=(pi/rlam)*sp*cp*xd+(1e-8);
yy=(pi/rlam)*sp*sp*yd+(1e-8);

fact5=(sin(xx)/xx)*(sin(yy)/yy);

ctcp=ct*cp;
ctsp=ct*sp;

hfact1=sp*hwk-cp*hbpol;
hfact3=ctcp*hwk+ctsp*hbpol;

htp=hfact1*hh*fact5*xd*yd*knum;
hpp=hfact3*hh*fact5*xd*yd*knum;

heth=zz0*hpp;
heph=-zz0*htp;

hex=ct*cp*heth-sp*heph;
hey=ct*sp*heth+cp*heph;
hez=-st*heth;




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