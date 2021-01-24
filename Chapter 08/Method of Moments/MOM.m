%C*****************************************************************
%C     MoM
%C*****************************************************************
%C     This is a MATLAB based program using
%C
%C       I.   POCKLINGTON'S INTEGRAL EQUATION (8-24)
%C       II.  HALLEN'S INTEGRAL EQUATION (8-27)
%C
%C     That computes:
%C
%C       A.  CURRENT DISTRIBUTION
%C       B.  INPUT IMPEDANCE 
%C       C.  NORMALIZED AMPLITUDE RADIATION PATTERN
%C
%C     Of a linear symmetrically-excited dipole. 
%C
%C     THIS PROGRAM USES PULSE EXPANSION FOR THE ELECTRIC CURRENT MODE 
%C     AND POINT-MATCHING FOR THE ELECTRIC FIELD AT THE CENTER OF EACH 
%C     WIRE SEGMENT.
%C
%C     DELTA-GAP FEED MODEL IS USED IN BOTH FORMULATIONS.  IN ADDITION,
%C     MAGNETIC-FRILL GENERATOR IS AVAILABLE IN THE POCKLINGTON'S
%C     INTEGRAL EQUATION.
%C
%C       OPTION I.   POCKLINGTON'S INTEGRAL EQUATION
%C       OPTION II.  HALLEN'S INTEGRAL EQUATION 
%C
%C                                                                     
%C         ** INPUT DATA:                                                       
%C
%C         TL  = TOTAL DIPOLE LENGTH (IN WAVELENGTHS)                 
%C         RA  = RADIUS OF THE WIRE (IN WAVELENGTHS)                        
%C         NM  = TOTAL NUMBER OF SUBSECTIONS (MUST BE AN ODD INTEGER)       
%C         IEX = OPTION TO USE EITHER MAGNETIC-FRILL GENERATOR OR DELTA GAP 
%C         
%C         IEX = 1 :  MAGNETIC-FRILL GENERATOR                            
%C         IEX = 2 :  DELTA-GAP FEED                                           
%C
%C         ** NOTE:  IGNORE INPUT PARAMETER IEX WHEN CHOOSING OPTION II
%C                   (i.e., HALLEN'S FORMULATION)
%C     *****************************************************************

clear;
close;
%*********************************************************
%***  OUTPUT DEVICE OPTION # : 1 ---> SCREEN
%                              2 ---> OUTPUT FILE
%*********************************************************
disp('OUTPUT DEVICE OPTION FOR THE OUTPUT PARAMETERS');
disp('        OPTION (1):SCREEN');
disp('        OPTION (2):OUTPUT FILE');
option_a=0; option_b=0; filename=0;
while((option_a~=1)&(option_a~=2))
   option_a=input('OUTPUT DEVICE =');
   if(option_a==2)
      filename=input('INPUT THE DESIRED OUTPUT FILENAME <in single quotes> = ','s');
   end
end
%*********************************************************
%***  READING OPTION # : 1 ---> POCKLINGTON'S EQUATIONS
%                        2 ---> HALLEN'S EQUATION
%*********************************************************
disp('CHOICE OF POCKLINGTON''S OR HALLEN''S EQN.');
disp('        OPTION (1):POCKLINGTON''S EQN.');
disp('        OPTION (2):HALLEN''S EQN.');
option=input('OPTION NUMBER =');
switch option
   
   
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
%***********************************************************************
%***           %POCKLINGTON'S INTEGRAL EQUATIONS                     ***
%***********************************************************************  
 


case 1
%******************************************************
%     SOME CONSTANTS AND ARRAYS 
%******************************************************
nmax=200;
nmt=2*nmax-1;
cgaw=zeros(nmax);
cga=cgaw(1,1:nmax);
zmnw=zeros(nmt);
zmn=zmnw(1,1:nmt);
waw=zeros(nmt);
wa=waw(1,1:nmt);
etm=zeros(361);
etmm=etm(1,1:361);
pi=3.14159265;
rad=pi/180;
beta=2.0*pi;
eta=120*pi;
%******************************************************
%***  INPUT DATA
%******************************************************
nm=input('NUMBER OF SUBDIVISIONS <ODD NUMBER> =');
tl=input('TOTAL DIPOLE LENGTH <WAVELENGTHS> =');
ra=input('RADIUS OF DIPOLE <WAVELENGTHS> =');
%*********************************************************
%***  EXITATION OPTION # : 1 ---> MAGNETIC-FRILL
%                          2 ---> DELTA-GAP
%**********************************************************
disp('EXCITATION :');
disp('        OPTION (1):MAGNETIC-FRILL');
disp('        OPTION (2):DELTA-GAP');
iex=input('OPTION NUMBER :');
hl=tl/2;
nmh=0.5*(nm+1);
dz=2*hl/nm;
zm=hl-0.5*dz;
b=0.5*dz;
a=-0.5*dz;
n=79;
hd=(b-a)/(n+1);
%*********************************************************************
%     THE IMPEDANCE MATRIX HAS A TOEPLITZ PROPERTY, THEREFORE ONLY
%     NM ELEMENTS NEED TO BE COMPUTED, AND THE MATRIX IS FILLED IN A
%     FORM THAT CAN BE SOLVED BY A TOEPLITZ MATRIX SOLVING 
%**********************************************************************
for I = 1: nm                
   		 zn=hl-(I-0.5)*dz;
   		 za1=zn-zm+a;
          %******************************************************************
          %     FAST ALGORITHM FORM OF THE SIMPSON'S INTEGRAL ROUTINE
          %******************************************************************
          recgp=sqrt(ra*ra+za1*za1);
          cgp1=exp(-j*beta*recgp)*((1.0+j*beta*recgp)*(2.0*recgp*recgp-3.0*ra*ra)+(beta*ra*recgp)^2)/(2.0*beta*recgp^5);                               
          zb1=zn-zm+b;
          roc=sqrt(ra*ra+zb1*zb1);
          cgp2=exp(-j*beta*roc)*((1.0+j*beta*roc)*(2.0*roc*roc-3.0*ra*ra)+(beta*ra*roc)^2)/(2.0*beta*roc^5);      
          crt=cgp1+cgp2;
          for k = 1: n
                xk=a+k*hd;
                zx1=zn-zm+xk;
                r=sqrt(ra*ra+zx1*zx1);
                cgp3=exp(-j*beta*r)*((1.0+j*beta*r)*(2.0*r*r-3.0*ra*ra)+(beta*ra*r)^2)/(2.0*beta*r^5);            
                if mod(k,2)~=0
                      crt=crt+4.0*cgp3;
                else
                      crt=crt+2.0*cgp3;
                end
            end   
        crt=crt*hd*0.33333;
        zmn(I)=crt;
      if I~=1
       zmn(nm+I-1)=crt;
      end
end
rb=2.3*ra;
tlab=2.0*log(2.3);
for i = 1: nm
      zi=hl-(i-0.5)*dz;
      r1=beta*sqrt(zi*zi+ra*ra);
      r2=beta*sqrt(zi*zi+rb*rb);
      j=0.0+j*1.0;
      if iex==1
         cga(i)=-j*beta^2/(eta*tlab)*(expm(-j*r1)/r1-expm(-j*r2)/r2);
      else  
            if i~=nmh
              cga(i)=0;
            else 
              cga(i)=-j*beta/(eta*dz);
            end
      end
   end
%******************************************************************
%       INPUT:
%       (C)A(2*M - 1)        THE FIRST ROW OF THE T-MATRIX FOLLOWED BY
%                            ITS FIRST COLUMN BEGINNING WITH THE SECOND
%                            ELEMENT.  ON RETURN A IS UNALTERED.
%       (C)B(M)              THE RIGHT HAND SIDE VECTOR B.
%       (C)WA(2*M-2)         A WORK AREA VECTOR
%       (I)M                 ORDER OF MATRIX A.
%     OUTPUT:
%       (C)B(M)              THE SOLUTION VECTOR.
%     PURPOSE:
%       SOLVE A SYSTEM OF EQUATIONS DESCRIBED BY A TOEPLITZ MATRIX.
%       A * X = B
%     ******************************************************************
t2=length(zmn);
t1=length(wa);
m=nm;
a=zmn;
b=cga;
wa=wa;
m=nm;
a1=a;
a2=a(m+1:t2);
t=b;
c1=wa;
c2=wa(m-1:t1);
r1=a1(1);
t(1)=b(1)/r1;
if m ~=1
    for n=2:m
        n1=n-1;
        n2=n-2;
        r5=a2(n1);
        r6=a1(n);
            if n ~= 2
                c1(n1)=r2;
                    for i1=1:n2
                        i2=n-i1;
                        r5=r5+a2(i1)*c1(i2);
                        r6=r6+a1(i1+1)*c2(i1);
                    end
            end
        r2=-r5/r1;
        r3=-r6/r1;
        r1=r1+r5*r3;
            if n~=2
                r6=c2(1);
                c2(n1)=0.0+j*0.0;
                    for i1=2:n1
                        r5=c2(i1);
                        c2(i1)=c1(i1)*r3+r6;
                        c1(i1)=c1(i1)+r6*r2;
                        r6=r5;
                    end
            end
        c2(1)=r3;
        r5=0.0+j*0.0;
            for i1=1:n1
                i2=n-i1;
                r5=r5+a2(i1)*t(i2);
            end
        r6=(b(n)-r5)/r1;
	        for i1=1:n1
                t(i1)=t(i1)+c2(i1)*r6;
            end
        t(n)=r6;
    end
 end
cga=t;
%******************************************************************
%     Output files for curr-MoM.dat
%******************************************************************
fid=fopen('Curr-MoM_m.dat','wt');
fprintf(fid,'CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE\n');
fprintf(fid,'POSITION Z:   MAGNITUDE:       REAL PART:       IMAGINARY PART:      SECTION #:\n\n');
%*****************************************************************************
%     OUTPUT THE CURRENT DISTRIBUTION ALONG OF THE DIPOLE
%*****************************************************************************
for i=1:nmh
    xi=hl-(i-0.5)*dz;
    yi=abs(cga(i));
    table2=[xi,yi,real(cga(i)),imag(cga(i)),i];
    fprintf(fid,'%1.4f       %1.6f            %1.6f        %1.6f          %2.1f \n\n',table2);
end
%******************************************************************
%     Output figure for curr-MoM.dat
%******************************************************************
i=[1:nmh];
xi=hl-(i-0.5)*dz;
yi=abs(cga(i));
stem(xi,yi)
grid on;
xlabel('Position (z)')
ylabel('Magnitude')
title('CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE')
figure;

fclose(fid);
%******************************************************
%     COMPUTATION OF THE INPUT IMPEDANCE
%******************************************************
zin=1.0/cga(nmh);
%******************************************************************
%     COMPUTATION OF AMPLITUDE RADIATION PATTERN OF THE ANTENNA
%******************************************************************
for i=1:181
    theta=(i-1.0)*rad;
    cth=cos(theta);
    sth=sin(theta);
      if abs(cth)<0.001
         ft=1;
      else
         ft=sin(beta*dz*cth*0.5)/(beta*dz*cth*0.5);
      end
    crt=0;
         for m=1:nm
             zm=hl-(m-0.5)*dz;
             crt=crt+exp(j*beta*zm*cth)*ft*cga(m)*dz;
         end
    ptt=abs(crt)*sth*sth*eta*0.5;
    etmm(i)=ptt;
 end
amax=etmm(1);
for i=2:181
    if etmm(i)>=amax
       amax=etmm(i);
    end
end
for i=1:181
    ptt=etmm(i)/amax;
      if ptt<=0.00001
         ptt=0.00001;
         etmm(i)=20*log10(ptt);
      end
      etmm(i)=20*log10(ptt);
end
%******************************************************************
%     Output files for Patt-MoM_m.dat
%******************************************************************
pok=fopen('Patt-MoM_m.dat','wt');
fprintf(pok,'RADIATION PATTERN   vs   OBSERVATION ANGLE THETA\n');
fprintf(pok,'THETA (in degrees)        MAGNITUDE(in dB)\n');

for i=1:181
  xi=i-1;
  table3=[xi,etmm(i)];
  fprintf(pok,'%3.1f                            %3.3f         \n',table3);
end
%******************************************************************
%     Output figure for Patt-MoM_m.dat
%******************************************************************
i=[1:181];
xi=i-1;
% Polar Plot
etmm1=[etmm(1:181),fliplr(etmm(1:180))];
q=polar_dB([0:360],etmm1,-40,0,4,'-');
set(q,'linewidth',1.5);
% plot(xi,etmm(i),'linewidth',2);
% axis ([0 180 -60 0]);
% grid on;
% xlabel('Theta(degrees)');
% ylabel('Magnitude(dB)');
title('RADIATION PATTERN   vs   OBSERVATION ANGLE')
for i=182:361
  xi=i-1;
   table4=[xi,etmm(362-i)];
   fprintf(pok,'%3.1f                           %3.3f         \n',table4);
end



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
%***********************************************************************
%***          HALLEN'S INTEGRAL EQUATION                             ***
%***********************************************************************

case 2
%******************************************************
%     SOME CONSTANTS AND INPUT DATA
%******************************************************
nm=input('NUMBER OF SUBDIVISIONS <ODD OR EVEN NUMBER> =');
tl=input('TOTAL DIPOLE LENGTH <WAVELENGTHS> =');
ra=input('RADIUS OF DIPOLE <WAVELENGTHS> =');
n=nm;
l=tl;
rho=ra;
rtod=180/pi;
eta=120*pi;
bk=2*pi;
cj=0.0+j*1.0;
dz=l/(2*(n-1));
dz1=l/(2*n);
%C***********************************************************************
%C    PERFORM COMPLEX INTEGRATION USING SIXTEEN POINT 
%C    GAUSSIAN QUADRATURE WITH INCREASING ACCURACY SET 
%C    BY INTEGER NO
%C***********************************************************************
absica=[-0.095012509837637,-0.281603550779259;
        -0.458016777657227,-0.617876244402644;
        -0.755404408355003,-0.865631202387832;
        -0.944575023073233,-0.989400934991650;
         0.095012509837637,0.281603550779259;
         0.458016777657227,0.617876244402644;
         0.755404408355003,0.865631202387832;
         0.944575023073233,0.989400934991650];
wght=[0.189450610455068,0.182603415044924;
         0.169156519395002,0.149595988816577;
         0.124628971255534,0.095158511682493;
         0.062253523938648,0.027152459411754;
         0.189450610455068,0.182603415044924;
         0.169156519395002,0.149595988816577;
         0.124628971255534,0.095158511682493;
         0.062253523938648,0.027152459411754];
      %************************************************************
      %     FILL THE MATRIX AND EXCITATION VECTOR OF THE SYSTEM: 
      %      [ZMATRX] and [ELECUR]
      %************************************************************
for i=1:n
    z=(2*i-1)*dz1/2.0;
    zmatrx(i,n)=-cos(bk*z);
    elecur(i)=-cj*sin(bk*z)/(2.0*eta);
    for j=1:n-1
        lower=(j-1)*dz1;
        upper=(j)*dz1;
        %*****************************************************************
        %     PERFORM NUMERICAL INTEGRATION OF THE KERNEL FOR HALLEN'S
        %     INTEGRAL EQUATION
        %*****************************************************************
        no=10;
        del=(upper-lower)/(2.0*no);
        sum=0.0+j*0.0;
        for t=1:no
            s=lower+(2*t-1)*del;
            for q=1:16
               x=s+absica(q)*del;
               %***********************************************************************
               %C    KERNEL PROVIDES THE KERNEL OF HALLEN'S EQN FOR INTEGRATION 
               %C    SYMMETRY IS USED TO REDUCE THE SYSTEM OF EQUATIONS AND 
               %C    HENCE ALTERS THE KERNEL
               %***********************************************************************
                r1=sqrt(rho*rho+(z-x)*(z-x));
                r2=sqrt(rho*rho+(z+x)*(z+x));
                kernel=exp(-cj*bk*r1)/(4.0*pi*r1)+exp(-cj*bk*r2)/(4.0*pi*r2);
                sum=sum+wght(q)*kernel;
            end
        end
        res=sum*del;
    zmatrx(i,j)=res;
  end
end

%*********************************************************************************
%   DECOMPOSE AND SOLVE THE SYSTEM FOR THE CURRENT DISTRIBUTION
%*********************************************************************************

%****************************
%    GET SCALING INFO.
%****************************
for i=1:n
    zmax=0.0;
    for j=1:n
        caz=abs(zmatrx(i,j));
        if caz >= zmax
            zmax=caz;
        end
    end
    scal(i)=1.0/zmax;
 end
%**************************** 
%    CROUT's algorithm.
%****************************
for j=1:n
    for i=1:j-1
        for k=1:i-1
            zmatrx(i,j)=zmatrx(i,j)-zmatrx(i,k)*zmatrx(k,j);
        end
    end
    zmax=0.0;
    %*****************************************
    %    SEARCH FOR LARGEST PIVOT ELEMENT.
    %*****************************************
    for i=j:n
        for k=1:j-1
            zmatrx(i,j)=zmatrx(i,j)-zmatrx(i,k)*zmatrx(k,j);
        end
        problem=scal(i)*abs(zmatrx(i,j));
        if problem >=zmax
            imax=i;
            zmax=problem;
        end
     end
     %***********************************
     %    INTERCHANGE THE ROWS.
     %***********************************
    if j~=imax
        for k=1:n
            temp=zmatrx(imax,k);
            zmatrx(imax,k)=zmatrx(j,k);
            zmatrx(j,k)=temp;
        end
        scal(imax)=scal(j);
    end
    iperm(j)=imax;
    %***********************************
    %    DIVIDE BY PIVOT ELEMENT.
    %***********************************
    if j~=n
        for i=j+1:n
            zmatrx(i,j)=zmatrx(i,j)/zmatrx(j,j);
        end
    end
 end
%****************************************************************** 
%    SOLVES LINEAR SYSTEM GIVEN THE LU DECOMPOSITION FROM LUDEC
%    FORCING VECTOR IS REPLACED WITH SOLUTION VECTOR UPON EXIT
%*****************************
%    FORWARD SUBSTITUTION.
%******************************************************************
for i=1:n
    temp=elecur(iperm(i));
    elecur(iperm(i))=elecur(i);
    for j=1:i-1
        temp=temp-zmatrx(i,j)*elecur(j);
    end
    elecur(i)=temp;
 end
%******************************** 
%    BACKWARD SUBSTITUTION.
%********************************
for i=1:n
    ii=n-i+1;
    temp=elecur(ii);
    for j=ii+1:n
        temp=temp-zmatrx(ii,j)*elecur(j);
    end
    elecur(ii)=temp/zmatrx(ii,ii);
 end
%*************************************************************************** 
%     COMPUTATION OF THE INPUT IMPEDANCE AND CURRENT DISTRIBUTION
%***************************************************************************
zin=1.0/elecur(1);
%******************************************************************
%     Output files for curr-MoM_m.dat
%******************************************************************
fid=fopen('Curr-MoM_m.dat','wt');
fprintf(fid,'CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE\n');
fprintf(fid,'POSITION Z:   CURRENT MAGNITUDE:    CURRENT PHASE:    SECTION:\n');
for i=1:n-1
   %******************************************************************
   %    THIS FUNCTION IS COMPUTES THE ARCTANGENT GIVEN X,Y.  IT IS 
   %    SIMILAR TO ATAN2 EXCEPT IT AVOIDS THE RUN TIME ERRORS ON 
   %    SOME MACHINES FOR SMALL ARGUMENTS.
   %******************************************************************
    smlt=0.0000001;
    if (abs(x)<=smlt) & (abs(y)<=smlt)
        btan2=0.0;
    else
        btan2=atan2(real(elecur(i)),imag(elecur(i)));
    end
    cur=abs(elecur(i));
    pha=rtod*btan2;
    table=[i*dz-dz/2.0,cur,pha,i];
    fprintf(fid,'%2.6f       %2.6f              %3.6f          %3.1f\n\n\n\n',table);
 end
%******************************************************************
%     Output figure for curr-MoM_m.dat
%******************************************************************

i=[1:n-1];
smlt=0.0000001;
    if (abs(x)<=smlt) & (abs(y)<=smlt)
        btan2=0.0;
    else
        btan2=atan2(real(elecur(i)),imag(elecur(i)));
    end
     cur=abs(elecur(i));
    pha=rtod*btan2;
stem(i*dz-dz/2.0,cur)
grid on
xlabel('Position (z)')
ylabel('Current magnitude')
title('CURRENT DISTRIBUTION ALONG ONE HALF OF THE DIPOLE')
fclose(fid);
figure;
%*********************************************************************************
%   CALCULATE THE RADIATION PATTERN OF THE ANTENNA
%*********************************************************************************
maxang=181;
pmax=-1.0;
for j=1:maxang
    theta=(j-1)/(maxang-1)*pi;
    %********************************************************************
    %    CALCULATES THE RADIATED POWER LEVEL AT ANGLE THETA RADIANS
    %    TO THE DIPOLE AXIS.  SINCE THE PATTERN IS NORMALIZED TO THE 
    %    MAXIMUM RADIATED POWER, COMMON CONSTANTS ARE REMOVED
    %********************************************************************
    sth=sin(theta);
    cth=cos(theta);
    arg=pi*dz*cth;
    if abs(arg) <= 0.001
        ft=1.0;
    else
        ft=sin(arg)/arg;
    end
    n2=n-1;;
    crt=0.0+j*0.0;
    for k=1:n2
        argp=pi*(-l+(k-1)*2.0*dz+dz)*cth;
        crt=crt+exp(cj*argp)*ft*elecur(k);
        argp=pi*(-l+(n2+k-1)*2.0*dz+dz)*cth;
        crt=crt+exp(cj*argp)*ft*elecur(k);
    end
    power=abs(crt)*sth*sth;
    pwr(j)=power;
    if pwr(j)>=pmax
        pmax=pwr(j);
    end
 end
 
%******************************************************************
%     Output files for Patt-MoM_m.dat
%******************************************************************
hel=fopen('Patt-MoM_m.dat','wt');
fprintf(hel,'RADIATION PATTERN   vs  OBSERVATION ANGLE THETA\n');
fprintf(hel,'THETA (in degrees)      MAGNITUDE (in dB)\n');
%*******************************************************************
%   WRITE THE RADIATION PATTERN IN dB
%*******************************************************************
for i=1:maxang
    theta=(i-1)/(maxang-1)*180;
    pattrn=pwr(i)/pmax;
    if pattrn <=0.00001
        pattrn=-100;
    else
        pattrn=20*log10(pattrn);
    end
    table1=[theta,pattrn];
    fprintf(hel,'%3.1f                         %3.3f         \n',table1);
end
%******************************************************************
%     Output figure for Patt-MoM_m.dat
%******************************************************************
i=[1:181];
theta=(i-1);
pattrn=pwr(i)/pmax;
pattrn=20*log10(pattrn);
% Polar Plot
pattrn=[pattrn,fliplr(pattrn(1:180))];
q=polar_dB([0:360],pattrn,-40,0,4,'-');
set(q,'linewidth',1.5);
% plot(theta,pattrn,'linewidth',2);
% axis ([0 180 -60 0]);
% grid on;
% xlabel('Theta(degrees)');
% ylabel('Magnitude(dB)');
title('RADIATION PATTERN   vs   OBSERVATION ANGLE')
for i=182:361
    theta=(i-1)/(maxang-1)*180;
    pattrn=pwr(362-i)/pmax;
    if pattrn <=0.00001
        pattrn=-100;
    else      
        pattrn=20*log10(pattrn);
    end
    table2=[theta,pattrn];
    fprintf(hel,'%3.1f                         %3.3f         \n',table2);
end
fclose(hel);
end
%*******************************************
%        OUTPUT FILE
%*******************************************
if(option_a==2)
   diary(filename);
end
%*******************************************
%         FORMAT STATEMENTS 
%*******************************************     
disp(sprintf('\n\n\n\n\n\n\n\n'));
disp(strvcat('WIRE ANTENNA PROBLEM','================'));
if (option==1)&(iex==1)
    disp(sprintf('POCKLINGTON''S EQUATION AND MAGNETIC FRILL MODEL'));
elseif (option==1)&(iex==2)
    disp(sprintf('POCKLINGTON''S EQUATION AND DELTA GAP MODEL'));
elseif (option==2)
    disp(sprintf('HALLEN''S EQUATION'));
end
disp(sprintf('\nLENGTH = %4.4f (WLS)',tl));
disp(sprintf('RADIUS OF THE WIRE = %4.4f (WLS)',ra));
disp(sprintf('NUMBER OF SUBSECTIONS = %2.2f\n',nm));
ZR=real(zin);
ZI=imag(zin);
if ZI>0
    output=sprintf('INPUT IMPEDANCE: Z= %5.1f +j %5.1f  (OHMS)\n',ZR,ZI);
else
    output=sprintf('INPUT IMPEDANCE: Z= %5.1f -j %5.1f  (OHMS)\n',ZR,abs(ZI));
end
disp(output);
disp(strvcat('*** NOTE:',...
   '    THE DIPOLE CURRENT DISTRIBUTION IS STORED IN Curr-MoM_m.dat',...
   '    THE AMPLITUDE RADIATION PATTERN IS STORED IN Patt-MoM_m.dat',...
   '    ========================================================='));
diary off;
warning off;