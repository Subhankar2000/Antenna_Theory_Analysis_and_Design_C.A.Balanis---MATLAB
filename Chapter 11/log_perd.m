function analdesg

%*****************************************
% written by                             %       
%                                        % 
% Tassos Panaretos                       %
% Arizona State University               %
% Tempe, AZ                              %
% 14 December 2002                       %
%*****************************************
% revised by                             % 
% Bo Yang                                %
% Arizona State University               %
% Tempe, AZ                              %
% 15 October 2004                        %
%*****************************************
clear all;
close all;

tic; 
   
%**************************************************************
%Initialize default settings
%**************************************************************

design_TITLE='Log-Periodic Dipole Array Design';  %1
% Fhigh       =216;                               %2
% Flow        =54;                               %3
% Tau         =0.865;                               %4
% sigm       =0.15825;                              %5
% D0          =8;                                %6
% LD          =146;                               %7
% Rs          =50.0;                                %8 
% LLin        =0.0;                                %9
% ZCin        =50.0+i*0.0;                         %10
% Rin         =0.0;                               %11
% LLout       =0.0;                                %12
% Zout        =50.0+i*0.0;                          %13
% AFhigh      =216;                               %14
% AFlow       =54;                                %15 
% AFpowr      =20;                                 %16
% AFSEH       =150;                               %17
% AFSC        =1500;                               %18
% Phi         =45;                                 %19 
% Navail      =2;                                 %20
% Davail( 1)  =0.47999;                             %21 
% Davail( 2)  =1.89999;                              %22
% Davail( 3)  =.1875;
% Davail( 4)  =.25;
% Davail( 5)  =.375;
% Davail( 6)  =.5;
% Davail( 7)  =.625;
% Davail( 8)  =.75;
% Davail( 9)  =.875;
% Davail(10)  =1.0;
% Davail(11)  =1.25;
% Davail(12)  =1.5;
% Davail(13)  =1.75;
% Davail(14)  =2.0;
% Davail(15)  =2.5;
% OutFlag     ='YYYNYYNO';                         %36  
% Quant       ='NY';                               %37 
% SB          =2.15899992;                         %38
% DB          =1.778;                              %39
    
%***********************************************

%------------------------------------------------------------------------------------------------  
% Design Frequency
%------------------------------------------------------------------------------------------------    
    ERR = 1;
    while(ERR ~= 0)
        Fhigh = str2num(input('Upper design frequency (in MHz) = ','s'));
        if (isempty(Fhigh))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Fhigh <= 0)
            fprintf('\n   *** ERROR: The number has to be greater than 0!\n\n');
        else
            ERR = 0;
        end
    end
    
    ERR = 1;
    while(ERR ~= 0)
        Flow = str2num(input('Lower design frequency (in MHz) = ','s'));
        if (isempty(Flow))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Flow <= 0)
            fprintf('\n   *** ERROR: The number has to be greater than 0!\n\n');
        elseif (Fhigh<Flow)
            fprintf('\n   *** ERROR: Lower frequency must be lower than higher frequency!\n\n');
        else
            ERR = 0;
        end
    end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Tau & Sigma / D0
%------------------------------------------------------------------------------------------------
figure;
[tausigma,cmap]=imread('tau_sigma.tif');
imagesc(tausigma); colormap(cmap); 
axis off;

ERR = 1;
while(ERR ~= 0)
    fprintf('\nSpecify either Tau and Sigma or the desired Directivity.\n');
    fprintf('Notice: Tau and Sigma are given by formulas (11-26) and (11-27) respectively.\n');
    fprintf('Above, there is the nomogram according to which, Tau and Sigma are specified for a given Directivity.\n');
    fprintf('If the Directivity is specified, it has to be greater than 7 dBi and smaller than 11 dBi.\n');
    fprintf('----------------------\n');
    fprintf('(1) Tau and Sigma\n');
    fprintf('(2) Desired directivity\n\n');
    SELECT_mode = str2num(input('SELECT = ','s'));
    if or((SELECT_mode == 1), (SELECT_mode == 2))
        ERR = 0;
    end
end

if SELECT_mode == 1
        ERR = 1;
    while(ERR ~= 0)
        Tau = str2num(input('Tau ( 0 <= Tau <= 1) = ','s'));
        if (isempty(Tau))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Tau < 0)|(Tau >1)
            fprintf('\n   *** ERROR: 0 <= Tau <= 1!\n\n');
        else
            ERR = 0;
        end
    end
    
        ERR = 1;
    while(ERR ~= 0)
        sigm = str2num(input('Sigma = ','s'));
        if (isempty(sigm))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (sigm < 0)
            fprintf('\n   *** ERROR: Sigma has to be greater than 0!\n\n');
        else
            ERR = 0;
        end
    end
    D0=0;
else
            ERR = 1;
    while(ERR ~= 0)
        D0 = str2num(input('Desired directivity (in dBi) = ','s'));
        if (isempty(D0))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (D0 < 7)|(D0 >11)
            fprintf('\n   *** ERROR: 7 <= D0 <= 11! For lower or higher directivity, please specify Tau and Sigma directly!\n\n');
        else
            ERR = 0;
        end
    end
[Tau,sigm]=FindST(D0);
end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Length-to-Diameter Ratio
%------------------------------------------------------------------------------------------------    
figure;
[drawing2,cmap]=imread('Drawing2.tif');
imagesc(drawing2); colormap(cmap); 
axis off;


            ERR = 1;
    while(ERR ~= 0)
        LD = str2num(input('Length-to-diameter ratio (LD > 20) = ','s'));
        if (isempty(LD))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (LD < 20)
            fprintf('\n   *** ERROR: LD > 20!\n\n');
        else
            ERR = 0;
        end
    end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Source resistance (in Ohms)
%------------------------------------------------------------------------------------------------    
            ERR = 1;
    while(ERR ~= 0)
        Rs = str2num(input('Source resistance (in Ohms) = ','s'));
        if (isempty(Rs))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Rs < 0)
            fprintf('\n   *** ERROR: Rs has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Length of source transmisssion line
%------------------------------------------------------------------------------------------------    
            ERR = 1;
    while(ERR ~= 0)
        LLin = str2num(input('Length of source transmisssion line (in m) = ','s'));
        if (isempty(LLin))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (LLin < 0)
            fprintf('\n   *** ERROR: LLin has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
                ERR = 1;
    while(ERR ~= 0)
        ZCin = str2num(input('Character impedance of source transmisssion line (in Ohms) = ','s'));
        if (isempty(ZCin))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (ZCin <= 0)
            fprintf('\n   *** ERROR: ZCin has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end

%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Length of termination transmisssion line
%------------------------------------------------------------------------------------------------    
            ERR = 1;
    while(ERR ~= 0)
        LLout = str2num(input('Length of termination transmisssion line (in m) = ','s'));
        if (isempty(LLout))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (LLout < 0)
            fprintf('\n   *** ERROR: LLout has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
    ERR=1;
    while(ERR ~= 0)
        Zout = str2num(input('Character impedance of termination transmisssion line (in Ohms) = ','s'));
        if (isempty(Zout))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Zout <= 0)
            fprintf('\n   *** ERROR: Zout has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end    
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Boom space and Boom diameter / Desired imput impedance and Boom diameter
%------------------------------------------------------------------------------------------------    
figure;
[drawing1,cmap]=imread('Drawing1.tif');
imagesc(drawing1); colormap(cmap); 
axis off;

ERR = 1;
while(ERR ~= 0)
    fprintf('\nSpecify boom spacing, boom diameter or desired input impedance, boom diameter.\n');
    fprintf('----------------------\n');
    fprintf('(1) Boom spacing and boom diameter\n');
    fprintf('(2) Desired input impedance, boom diameter\n\n');
    SELECT_mode = str2num(input('SELECT = ','s'));
    if or((SELECT_mode == 1), (SELECT_mode == 2))
        ERR = 0;
    end
end

if SELECT_mode == 1
                ERR = 1;
    fprintf('Note: Boom spacing must be greater than boom diameter.\n');
    while(ERR ~= 0)
        SB = str2num(input('Boom spacing(in cm) = ','s'));
        if (isempty(SB))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (SB <= 0)
            fprintf('\n   *** ERROR: SB has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
                    ERR = 1;
    while(ERR ~= 0)
        DB = str2num(input('Boom diameter(in cm) = ','s'));
        if (isempty(DB))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (DB <= 0)
            fprintf('\n   *** ERROR: Boom diameter has to be greater than zero!\n\n');
        elseif (DB > SB)
            fprintf('\n   *** ERROR: Boom diameter has to be smaller than boom spacing!\n\n');
        else
            ERR = 0;
        end
    end
    Rin=0;
else
            ERR = 1;
    while(ERR ~= 0)
        Rin = str2num(input('Desired input resistance(in Ohms) = ','s'));
        if (isempty(Rin))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Rin <= 0)
            fprintf('\n   *** ERROR: Rin has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
                    ERR = 1;
    while(ERR ~= 0)
        DB = str2num(input('Boom diameter(in cm) = ','s'));
        if (isempty(DB))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (DB <= 0)
            fprintf('\n   *** ERROR: Boom diameter has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
    SB=0;    
end

%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Tube diameter choice
%------------------------------------------------------------------------------------------------    
ERR = 1;
while(ERR ~= 0)
    fprintf('Round element diameters to the nearest size?\n');
    fprintf('----------------------\n');
    fprintf('(1) Yes\n');
    fprintf('(2) No\n\n');
    SELECT_mode = str2num(input('SELECT = ','s'));
    if or((SELECT_mode == 1), (SELECT_mode == 2))
        ERR = 0;
    end
end

if SELECT_mode ==1
    Quant(2)='Y';
    ERR=1;
    while(ERR ~= 0)
        Navail = str2num(input('Number of available tube diameters = ','s'));
        if (isempty(Navail))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Navail <= 0)
            fprintf('\n   *** ERROR: Number has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end
    
    fprintf('Note: Please enter the diameter from the smaller one to the larger one.\n');

    temp1=0;
    for i=1:Navail
        ERR=1;
            while(ERR ~= 0)
        temp = str2num(input(['Diameter No.',num2str(i),'(in cm) = '],'s'));
        if (isempty(temp))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (temp <= 0)
            fprintf('\n   *** ERROR: Number has to be greater than zero!\n\n');
        elseif (temp <= temp1)
            fprintf('\n   *** ERROR: Current diameter has to be greater than previous diameter!\n\n');
        else
            ERR = 0;
            Davail(i)=temp;
        end
        temp1=temp;        
    end
    end
else
    Quant(2)='N';
end    
%------------------------------------------------------------------------------------------------    

% Outflag
OutFlag='YYYYYYYY';

%------------------------------------------------------------------------------------------------    
% Calculate E- and H-plane patterns
%------------------------------------------------------------------------------------------------    
            ERR = 1;
    while(ERR ~= 0)
        AFSEH = str2num(input('Analysis frequency to calculate E- and H- plane patterns (in MHz) = ','s'));
        if (isempty(AFSEH))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (AFSEH < Flow)|(AFSEH >Fhigh)
            fprintf('\n   *** ERROR: Analysis frequency has to be between lower and upper frequency!\n\n');
        else
            ERR = 0;
        end
    end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Custom pattern plot
%------------------------------------------------------------------------------------------------    
            ERR = 1;
    while(ERR ~= 0)
        Phi = str2num(input('Phi to calculate custom patterns (in degrees) = ','s'));
        if (isempty(Phi))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (Phi < 0)|(Phi > 360)
            fprintf('\n   *** ERROR: 0 <= Phi <= 360 degrees!\n\n');
        else
            ERR = 0;
        end
    end

            ERR = 1;
    while(ERR ~= 0)
        AFSC = str2num(input('Analysis frequency to calculate custom patterns (in MHz) = ','s'));
        if (isempty(AFSC))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (AFSC < Flow)|(AFSC >Fhigh)
            fprintf('\n   *** ERROR: Analysis frequency has to be between lower and upper frequency!\n\n');
        else
            ERR = 0;
        end
    end
%------------------------------------------------------------------------------------------------    

%------------------------------------------------------------------------------------------------    
% Gain vs. frequency
%------------------------------------------------------------------------------------------------  
    fprintf('Swept frequency analysis:\n');
    fprintf('--------------------------\n');
            ERR = 1;
    while(ERR ~= 0)
        AFhigh = str2num(input('Upper analysis frequency (in MHz) = ','s'));
        if (isempty(AFhigh))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (AFhigh > Fhigh)
            fprintf('\n   *** ERROR: Analysis frequency has to be less than Fhigh!\n\n');
        else
            ERR = 0;
        end
    end

                ERR = 1;
    while(ERR ~= 0)
        AFlow = str2num(input('Lower analysis frequency  (in MHz) = ','s'));
        if (isempty(AFlow))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (AFlow < Flow)|(AFlow > AFhigh)
            fprintf('\n   *** ERROR: Flow <= AFlow <= AFhigh <= Fhigh!\n\n');
        else
            ERR = 0;
        end
    end

                ERR = 1;
    while(ERR ~= 0)
        AFpowr = str2num(input('Number of frequency steps per octave = ','s'));
        if (isempty(AFpowr))
            fprintf('\n   *** ERROR: Please enter a number!\n\n');
        elseif (AFpowr <= 0)
            fprintf('\n   *** ERROR: AFpowr has to be greater than zero!\n\n');
        else
            ERR = 0;
        end
    end

%------------------------------------------------------------------------------------------------    


% Main program
%*************************************************
ms = 30;             % max number of elements
C  = 299.7925;       % speed of light in in Mkm/sec

%*********************************************

Alpha=atan((1-Tau)/4/sigm);

Bs=(Fhigh/Flow)*(1.1+(7.7*(1-Tau)^2)/tan(Alpha));

Lmax=C/Flow;

N=floor( 1-log(Bs)/log(Tau) );

L(N)=Lmax/2;
ZL(N)=L(N)/2/tan(Alpha);
D(N)=L(N)/LD/.01;

for ii=N-1:-1:1
  L(ii)=L(ii+1)*Tau;
  ZL(ii)=ZL(ii+1)*Tau;
  D(ii)=D(ii+1)*Tau;
end

%*****************************************************************
% QUANTIZE ELEMENTS TO NEAREST TUBE DIAMETER
% AND FIND TRUE L/D RATIO, Kavg
%*****************************************************************
  Dsave=0;
if upper(Quant(2))=='Y'
  Kavg=0;
  for jj=1:N,
    if (Navail~=0)
      Dbar=1.e6;
      for ii=1:Navail,
        if abs(D(jj)-Davail(ii))<Dbar
          Dbar=abs(D(jj)-Davail(ii));
          Dsave=Davail(ii);
%          keyboard
        end
      end
%      keyboard
      D(jj)=Dsave;
      Kavg=Kavg+L(jj)/D(jj)/N/.0254;
    else
      Kavg=LD;
    end
  end
else
  Kavg=LD;
end

%******************************************************************
% STEPS 5. - 7.  Calculate SB if it is not specified
%******************************************************************

Z0=0;

if SB==0
  
  ZA=120*(log(Kavg)-2.25);

  Sigmap=sigm/sqrt(Tau);
  dummy1=1/(8*Sigmap*ZA/Rin);

  Z0=Rin*(dummy1+sqrt(dummy1+1));
 
  SB=DB*cosh(Z0/120);

end

%*****************************************************************
%
% MAKE A DESIGN SUMMARY IF REQUIRED
%
%*****************************************************************

disp('Writting file LP_DES.OUT  -- Design Summary')

if upper(OutFlag(1))=='Y'
  
   output2(N, ZL, L, D, Z0,design_TITLE, Fhigh, Flow,...
                 Tau, sigm, D0, LD, Rs, LLin,ZCin,...
                  Rin, LLout, Zout, SB, DB)

end

%****************************************************************
%
% ANALYZE THE ANTENNA
%
%****************************************************************
 
% if upper(OutFlag(5))=='Y'
  
  if exist('LP_swept.dat')~=0
     disp('Deleting file LP_swept.dat')
     if (isunix)
        ! rm LP_swept.dat
     else
        ! del LP_swept.dat
     end   
  end

  if exist('LP_gain.dat')~=0
    disp('Deleting file LP_gain.dat')
    if (isunix)
       ! rm LP_gain.dat
    else
       ! del LP_gain.dat
    end   
  end

  if exist('LP_ftbr.dat')~=0
     disp('Deleting file LP_ftbr.dat')
     if (isunix)
        ! rm LP_ftbr.dat
     else
        ! del LP_ftbr.dat
     end   
  end

  if exist('LP_z_src.dat')~=0
     disp('Deleting file LP_z_src.dat')
     if (isunix)
        ! rm LP_z_src.dat
     else
        ! del LP_z_src.dat
     end
  end

  if exist('LP_z_ant.dat')~=0
     disp('Deleting file LP_z_ant.dat')
     if (isunix)
        ! rm LP_z_ant.dat
     else
        ! del LP_z_ant.dat
     end   
  end

  if exist('LP_vswr.dat')~=0
     disp('Deleting file LP_vswr.dat')
     if (isunix)
        ! rm LP_vswr.dat
     else
        ! del LP_vswr.dat
     end  
  end

% end

WATTS=1;

%****************************************************************
%
% DO CUSTOM PLANE ANALYSIS
%
%****************************************************************

if (upper(OutFlag(4))=='Y')
  OutFlag(8)='4';
  FMHZ=AFSC;
  disp('Custom plane analysis')
  disp(['Analysis Frequency in MHz ' num2str(FMHZ)])

  analyze(design_TITLE,N,ZL,LLin,LLout,L,D,DB,SB,...
	     Zout,ZCin,Rs,FMHZ,OutFlag,WATTS,Fhigh,Flow,Tau,sigm,Phi);

end



%****************************************************************
%
% DO E PLANE ANALYSIS
%
%****************************************************************

if (upper(OutFlag(2))=='Y')
  OutFlag(8)='2';
  FMHZ=AFSEH;
  disp('E plane analysis     ')
  disp(['Analysis Frequency in MHz ' num2str(FMHZ)])

  analyze(design_TITLE,N,ZL,LLin,LLout,L,D,DB,SB,...
          Zout,ZCin,Rs,FMHZ,OutFlag,WATTS,Fhigh,Flow,Tau,sigm,Phi)
end


%****************************************************************
%
% DO H PLANE ANALYSIS
%
%****************************************************************

if (upper(OutFlag(3))=='Y')
  OutFlag(8)='3';
  FMHZ=AFSEH;
  disp('H plane analysis     ')
  disp(['Analysis Frequency in MHz ' num2str(FMHZ)])

  analyze(design_TITLE,N,ZL,LLin,LLout,L,D,DB,SB,...
          Zout,ZCin,Rs,FMHZ,OutFlag,WATTS,Fhigh,Flow,Tau,sigm,Phi)

end

% SWEPT FREQUENCY ANALYSIS 

if OutFlag(5)=='Y'
  OutFlag(8)='5';
  AFstep=2^(1/AFpowr);
  FMHZ=AFlow;

  disp('Writting file LP_SWEPT.DAT -- Composite file ');
  disp('Writting file LP_GAIN.DAT  -- Gain vs. Frequency');
  disp('Writting file LP_FTBR.DAT  -- FRONT TO BACK RATIO');
  disp('Writting file LP_Z_SRC.DAT -- Impedance at Source');
  disp('Writting file LP_Z_ANT.DAT -- Impedance at Antenna');
  disp('Writting file LP_VSWR.DAT  -- VSWR in Source Line');

  fid_222=fopen('LP_swept.dat','a');

  fprintf(fid_222,'      Freq           Peak Dir.   Gain on Boresight  F/B        F/SLL        Re(ZinA)       Im(ZinA)       Mag(ZinA)     Re(ZinAS)     Im(ZinAS)     Mag(ZinAS)      VSWR');
  fprintf(fid_222,'\n');  
  fprintf(fid_222,'       MHz             dBi             dbi           db          db           ohms           ohms           ohms          ohms          ohms           ohms         ');

 disp('Swept frequency analysis')

% DO ANALYSIS
h=waitbar(0,'Program is running ...');
AFrange=AFhigh-AFlow;
  while FMHZ<AFhigh
    disp(['Frequency (MHz) : ' num2str(FMHZ)]);
    analyze(design_TITLE,N,ZL,LLin,LLout,L,D,DB,SB,Zout,ZCin,Rs,FMHZ,OutFlag,WATTS,Fhigh,Flow,...
            Tau,sigm,Phi);
    FMHZ=FMHZ*AFstep;
    waitbar((FMHZ-AFlow)/AFrange,h);
  end
end   
close(h);

%*********************************************************************
if OutFlag(2)=='Y'

 temp=load('LP_epat.dat');
 Epat=temp(:,2);
 Epat=Epat-max(Epat);
 
 temp=load('LP_hpat.dat'); 
 Hpat=temp(:,2);
 Hpat=Hpat-max(Hpat);
 
 theta=[1:360]';
 theta=theta-180;
 figure;
 polar_dB(theta,Epat,-60,0,6,'-');
 title('Normalized E-plane pattern'); 

 figure;
 polar_dB(theta,Hpat,-60,0,6,'-');
 title('Normalized H-plane pattern'); 

end

if OutFlag(4)=='Y'

 temp=load('LP_cpat.dat'); 
 Cpat=temp(:,2);
 Cpat=Cpat-max(Cpat);

 theta=[1:360]';
 theta=theta-180;
 figure;
 polar_dB(theta,Cpat,-60,0,6,'-');
 title(['Normalized pattern at custom-plane(\phi = ',num2str(Phi),' degrees)']); 

end

if OutFlag(5)=='Y'

 LPFTBR=load('LP_ftbr.dat');
 LPGAIN=load('LP_gain.dat');
%  LPSWEPT=load('LP_swept.dat');
 LPVSWR=load('LP_vswr.dat');
 LPZANT=load('LP_z_ant.dat');
 LPZSRC=load('LP_z_src.dat');

 figure;
 plot(LPFTBR(:,1),LPFTBR(:,2));
 xlabel('Frequency (MHz)');
 ylabel('Front-to-Back Ratio (dB)');
 grid on;
 xlim([AFlow AFhigh]);

 figure;
 plot(LPGAIN(:,1),LPGAIN(:,2));
 xlabel('Frequency (MHz)');
 ylabel('Gain (dBi)');
 grid on;
 xlim([AFlow AFhigh]);

 figure;
 plot(LPVSWR(:,1),LPVSWR(:,2));
 xlabel('Frequency (MHz)');
 ylabel('VSWR');
 grid on;
 xlim([AFlow AFhigh]);

 figure;
 plot(LPZANT(:,1),LPZANT(:,2));
 xlabel('Frequency (MHz)');
 ylabel('Input Impedance at Antenna (Ohms)');
 grid on;
 xlim([AFlow AFhigh]);
 
 figure;
 plot(LPZSRC(:,1),LPZSRC(:,2));
 xlabel('Frequency (MHz)');
 ylabel('Input Impedance at Source (Ohms)');
 grid on;
 xlim([AFlow AFhigh]);
  
end
%**********************************************************************
fclose all;
toc



%**********************************************************************
%**********************************************************************
%**********************************************************************

function analyze(design_TITLE, N, ZL, LLin, LLout, L, D, DB, SB,...
                    Zout, ZCin, Rs, FMHZ, OutFlag, WATTS, Fhigh,...
                    Flow, Tau, sigm, Phi)

NUMPOINT = 360;
big      = 1.e12;
c        = 299.7925;

BETA=2*pi*FMHZ/c;

Z0=120*log( SB/DB + sqrt( (SB/DB)^2-1 ) );

D=D*.01;

%disp('OK 1')
%keyboard
ZA=find_ZA(N,BETA,D,ZL,L);       

%disp('OK 2')
%keyboard

D=D/.01;

DZL(N)=LLout;
DZL(1:N-1)=ZL(2:N)-ZL(1:N-1);
DZL_0=LLin;

% EXCITATION MATRIX 

I=[1; zeros(N-1,1)];

YT=find_YT(N,DZL,BETA,Z0,Zout);

%YA=inv(ZA);

[L1 U1 P1]=lu(ZA);

YA=inv(U1)*inv(L1)*P1;

YAYT=YA+YT;

[L2 U2 P2]=lu(YAYT);

ZAZT=inv(U2)*inv(L2)*P2;

%ZAZT=inv(YAYT);

Vel=ZAZT*I;

Iel=YA*Vel;

% TERMINATION CURRENT

ZR=Z0 * (Zout +i * Z0   * tan(BETA*LLout) )/...
        (Z0   +i * Zout * tan(BETA*LLout) );

if abs(ZR)~=0
  Iout=Vel(N)/ZR;
else
  Iout=big;
end

Vout = Vel(N)*cos(BETA*LLout)-i*Iout*Z0*sin(BETA*LLout);
Iout = Iout * cos(BETA*LLout)-i*Vel(N)/Z0*sin(BETA*LLout);

ZinA=Vel(1)/1;

% VSWR

ABS_GAM=abs( (ZinA-ZCin)/(ZinA+ZCin) );

if ABS_GAM~=1
  VSWR=(1+ABS_GAM)/(1-ABS_GAM);
else
  VSWR=100000;
end

Iin = cos(BETA*LLin)+i*Vel(1)*sin(BETA*LLin)/ZCin;
Vin = Vel(1)*cos(BETA*LLin)+i*ZCin*sin(BETA*LLin);
ZinAS=Vin/Iin;

Vin     = (Rs+ZinAS)/ZinAS*Vin;

POWERIN = abs(Iin)^2*(real(ZinAS)+Rs)/2;

SCALE   = sqrt(WATTS/POWERIN);

Iel     = Iel*SCALE;
Vel     = Vel*SCALE;

Iin     = Iin*SCALE;
Vin     = Vin*SCALE;
Iout    = Iout*SCALE;
Vout    = Vout*SCALE;

RADIOUS = 100/BETA*2*pi;

Gain_Iso=10*log10(WATTS/4/pi/RADIOUS^2);

% CUSTOM PLANE CUT

if OutFlag(8)=='4'

  for kk=1:NUMPOINT,
    Theta     = kk*360/NUMPOINT;
    CPAT(kk)  = pattern(BETA,Phi,Theta,L,ZL,N,Iel,RADIOUS);
    ANGLE(kk) = Theta;
  end

% FIND FTBR FSLL

  [PeakSLL MainPeak BackLobe]=Find_SLL(NUMPOINT,CPAT);

  FSLL = CPAT(NUMPOINT/2) - PeakSLL;
  FTBR = CPAT(NUMPOINT/2) - BackLobe;
  PD0  = MainPeak         - Gain_Iso;
  D0A  = CPAT(NUMPOINT/2) - Gain_Iso;

end

% Make E-plane cut 

if (OutFlag(8)=='2')|(OutFlag(8)=='5')

  Phi=90;
  for kk=1:NUMPOINT,
    Theta     = kk*360/NUMPOINT;
    EPAT(kk)  = pattern(BETA,Phi,Theta,L,ZL,N,Iel,RADIOUS);
    ANGLE(kk) = Theta;
  end

% FIND FTBR FSLL

  [PeakSLL MainPeak BackLobe]=Find_SLL(NUMPOINT,EPAT);

  FSLL = EPAT(NUMPOINT/2) - PeakSLL;
  FTBR = EPAT(NUMPOINT/2) - BackLobe;
  PD0  = MainPeak         - Gain_Iso;
  D0A  = EPAT(NUMPOINT/2) - Gain_Iso;
  
  if exist('out_var.mat')&(isunix)
     ! rm out_var.mat 
     save out_var.mat D0A FTBR ZinAS VSWR
  
  elseif exist('out_var.mat')&(~isunix)
     ! del out_var.mat 
     save out_var.mat D0A FTBR ZinAS VSWR
     
  elseif (~exist('out_var.mat'))
     save out_var.mat D0A FTBR ZinAS VSWR
     
  end
  
  
end

% MAKE H-PLANE CUT

if OutFlag(8)=='3'
    
  Phi=0;
  for kk=1:NUMPOINT,  
    Theta     = kk*360/NUMPOINT;
    HPAT(kk)  = pattern(BETA,Phi,Theta,L,ZL,N,Iel,RADIOUS);
    ANGLE(kk) = Theta;
  end

% FIND FTBR FSLL

  [PeakSLL MainPeak BackLobe]=Find_SLL(NUMPOINT,HPAT);

  FSLL = HPAT(NUMPOINT/2) - PeakSLL;
  FTBR = HPAT(NUMPOINT/2) - BackLobe;
  PD0  = MainPeak         - Gain_Iso;
  D0A  = HPAT(NUMPOINT/2) - Gain_Iso;

end

Vin  = abs(Vin) + i*angle(Vin)*180/pi;
Iin  = abs(Iin) + i*angle(Iin)*180/pi;
Vout = abs(Vout)+ i*angle(Vout)*180/pi;
Iout = abs(Iout)+ i*angle(Iout)*180/pi;

if real(Vin)~=0
  Vin = 20*log10(real(Vin)) + i*imag(Vin);
else
  Vin = -1000 + i*mag(Vin);
end

if real(Iin)~=0
  Iin = 20*log10(real(Iin)) + i*imag(Iin);
else
  Iin = -1000 + i*mag(Iin);
end

if real(Vout)~=0
  Vout = 20*log10(real(Vout)) + i*imag(Vout);
else
  Vout =-1000 + i*mag(Vout);
end
  
if real(Iout)~=0
  Iout = 20*log10(real(Iout)) + i*imag(Iout);
else
  Iout =-1000 + i*mag(Iout);
end

for kk=1:N,

  Vel(kk) = abs(Vel(kk))+ i*angle(Vel(kk))*180/pi;
  Iel(kk) = abs(Iel(kk))+ i*angle(Iel(kk))*180/pi;

  if real(Vel(kk))~=0
    Vel(kk) = 20*log10(real(Vel(kk))) + i*imag(Vel(kk));
  else
    Vel(kk) = -1000 + i*imag(Vel(k));
  end

  if real(Iel(kk))~=0
    Iel(kk) = 20*log10(real(Iel(kk)))+i*imag(Iel(kk));
  else
    Iel(kk) = -1000 + i*imag(Iel(kk));
  end
  
end

% WRITE SWEPT FREQUENCY DATA

if OutFlag(8)=='5'

  fid_2 =fopen('LP_swept.dat','a'); 
  fprintf(fid_2,'\n %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f  %12.5f',FMHZ,PD0,D0A,FTBR,FSLL,real(ZinA),imag(ZinA),abs(ZinA),real(ZinAS),imag(ZinAS),abs(ZinAS),VSWR);
  fclose(fid_2);

  fid_2 =fopen('LP_gain.dat','a'); 
  fprintf(fid_2,'\n %3.5f  %3.5f',FMHZ,D0A);
  fclose(fid_2);

  fid_2 =fopen('LP_ftbr.dat','a');
  fprintf(fid_2,'\n %3.5f  %3.5f',FMHZ,FTBR);
  fclose(fid_2);

  fid_2 =fopen('LP_z_ant.dat','a');
  fprintf(fid_2,'\n %3.5f  %3.5f',FMHZ,abs(ZinA));
  fclose(fid_2);

  fid_2 =fopen('LP_z_src.dat','a');
  fprintf(fid_2,'\n %3.5f  %3.5f',FMHZ,abs(ZinAS));
  fclose(fid_2);
  
  fid_2 =fopen('LP_vswr.dat','a');
  fprintf(fid_2,'\n %3.5f  %3.5f',FMHZ,VSWR);
  fclose(fid_2);
  
end

%WRITE SINGLE FREQUENCY DATA

if OutFlag(8)~='5'
  
  if OutFlag(8)=='2'
    disp('Writting file LP_EPL.OUT -- Analysis in E- plane')
  end
 
  if OutFlag(8)=='3'
    disp('Writting file LP_HPL.OUT -- Analysis in H- plane')
  end

  if OutFlag(8)=='4'
    disp('Writting file LP_CPL.OUT -- Analysis in custom plane')
  end

  output1(N, Rs, ZL, LLin, LLout, L, D, DB, SB, Zout,...
         FMHZ, Vin, Iin, Vout, Iout, WATTS, ZCin, Iel, Vel,...
         ZinA, ZinAS, Z0, VSWR, PD0, D0A, FTBR, FSLL,design_TITLE,... 
         OutFlag, Fhigh, Flow, Tau, sigm, Phi);


  if OutFlag(8)=='2'

    if exist('LP_epat.dat')~=0
       disp('Deleting file LP_epat.dat')
       if (isunix)
          ! rm LP_epat.dat
       else
          ! del LP_epat.dat
       end  
    end

    fid_33=fopen('LP_epat.dat','wt');
    fprintf(fid_33,'%% Angle (degrees)   E-plane Gain (dB)');
    
    disp('Writting file LP_EPAT.DAT - E plane gain')
    for kk=1:NUMPOINT,
      Theta=kk*360/NUMPOINT;
      fprintf(fid_33,'\n%7.2f              %10.3f',Theta,EPAT(kk)-Gain_Iso);
    end
  
    fclose(fid_33);

  end


  if OutFlag(8)=='3'

    if exist('LP_hpat.dat')~=0
       disp('Deleting file LP_hpat.dat')
       if (isunix)
          ! rm LP_hpat.dat
       else
          ! del LP_hpat.dat
       end   
    end

    fid_33=fopen('LP_hpat.dat','wt');
    fprintf(fid_33,'%% Angle (degrees)   H-plane Gain (dB)');
    
    disp('Writting file LP_HPAT.DAT - H plane gain')
    for kk=1:NUMPOINT,
      Theta=kk*360/NUMPOINT;
      fprintf(fid_33,'\n%7.2f              %10.3f',Theta,HPAT(kk)-Gain_Iso);
    end
  
    fclose(fid_33);

  end

  if OutFlag(8)=='4'

    if exist('LP_cpat.dat')~=0
       disp('Deleting file LP_cpat.dat')
       if (isunix)
          ! rm LP_cpat.dat
       else
          ! del LP_cpat.dat
       end   
    end

    fid_33=fopen('LP_cpat.dat','wt');
    fprintf(fid_33,'%% Angle (degrees)   Custom-plane Gain (dB)');

    disp('Writting file LP_CPAT.DAT - H plane gain')
    for kk=1:NUMPOINT,
      Theta=kk*360/NUMPOINT;
       fprintf(fid_33,'\n%7.2f              %10.3f',Theta,CPAT(kk)-Gain_Iso);
      %fprintf(fid_33,'\n%3.2f %5.3f',Theta,CPAT(kk)-Gain_Iso);
    end
  
    fclose(fid_33);

  end
end

%  
%******************************************************************
%******************************************************************

function [relmax,absmax,backlobe]=Find_SLL(N,Gain);

 Nover2=N/2;
 backlobe=-10000;
 absmax  =-10000;
 relmax  =-10000;

for kk=1:N,

   peak=0;
   if ( ( Gain(mod(kk-1+N,N)+1) > Gain(mod(kk-2+N,N)+1))&... 
        ( Gain(mod(kk-1+N,N)+1) >=Gain(mod(kk    ,N)+1))&...
        ( Gain(mod(kk-1+N,N)+1) > Gain(mod(kk+1  ,N)+1)) )
      peak=1;
   end
     
   if ( (Gain(mod(kk-1+N,N)+1)> Gain(mod(kk-2+N,N)+1))&... 
        (Gain(mod(kk-1+N,N)+1)> Gain(mod(kk    ,N)+1)))
      peak=1;
   end

   if peak==1
     if (Gain(mod(kk-1,N)+1)>=absmax)
       relmax  = absmax;
       absmax  = Gain(kk);
       backlobe= Gain(mod(kk+Nover2-1,N)+1);
     elseif Gain(mod(kk-1,N)+1) > relmax
       relmax  = Gain(mod(kk-1,N)+1);
     end
   end 

end


%******************************************************************
%******************************************************************

function [Tau,sigm]=FindST(D0)

%******************************************************************
% This routine finds Sigma and Tau given the desired directivity.
% This table represents a linear interpolation of the points on the
% optimum directivity line of Figure 11.13. CAB Antennas: 
% Analysis and Design, edition 2
%******************************************************************

  if ((D0>= 7.0)&(D0< 7.5))     
     Tau=0.78 +(D0- 7.0)/0.5*(0.044);
  elseif ((D0>= 7.5)&(D0< 8.0))   
     Tau=0.824+(D0- 7.5)/0.5*(0.041);
  elseif ((D0>= 8.0)&(D0< 8.5)) 
     Tau=0.865+(D0- 8.0)/0.5*(0.030);
  elseif ((D0>=8.5)&(D0< 9.0))   
     Tau=0.895+(D0- 8.5)/0.5*(0.022);
  elseif ((D0>=9.0)&(D0<9.5)) 
     Tau=0.917+(D0- 9.0)/0.5*(0.008);
  elseif ((D0>=9.5)&(D0<10.0))      
     Tau=0.925+(D0- 9.5)/0.5*(0.017);
  elseif ((D0>=10.0)&(D0<10.5)) 
     Tau=0.942+(D0-10.0)/0.5*(0.013);
  elseif ((D0>10.5)&(D0<11.0)) 
     Tau=0.955+(D0-10.5)/0.5*(0.012);
end

sigm = 0.237838*Tau - 0.047484;
      

%***************************************************************
%***************************************************************

function output1(N,Rs,ZL,LLin,LLout,L,D,DB,SB,Zout,FMHZ,Vin,Iin,Vout,Iout,...
                 WATTS,ZCin,IA,VA,ZinA,ZinAS,Z0,VSWR,PD0,D0A,FTBR,FSLL,...
                 design_TITLE,OutFlag,Fhigh,Flow,Tau,sigm,Phi)


if OutFlag(8)=='2'
  fid_22=fopen('LP_Epl.out','wt');
end

if OutFlag(8)=='3'
  fid_22=fopen('LP_Hpl.out','wt');
end

if OutFlag(8)=='4'
  fid_22=fopen('LP_Cpl.out','wt');
end

fprintf(fid_22,'DIPOLE ARRAY DESIGN');
fprintf(fid_22,'\n');
fprintf(fid_22,design_TITLE);
fprintf(fid_22,'\n');
fprintf(fid_22,'\n');

if OutFlag(8)=='2'
  fprintf(fid_22,'E-plane Pattern Data');
  fprintf(fid_22,'\n');
end

if OutFlag(8)=='3'
  fprintf(fid_22,'H-plane Pattern Data');
  fprintf(fid_22,'\n');
end

if OutFlag(8)=='4'
  fprintf(fid_22,'Custom Plane Pattern Data');
  fprintf(fid_22,'\n');
  fprintf(fid_22,'Plane Angle =  %3.2f degrees',Phi);
end

fprintf(fid_22,'\n');
fprintf(fid_22,'Ele.        Z        L       D      Volts              Current           ');
fprintf(fid_22,'\n');
fprintf(fid_22,'           (m)     (m)     (cm)     (dBV)     (deg)    (dBAmp)     (deg)');
fprintf(fid_22,'\n');  

R1=real(Vin);
I1=imag(Vin);
R2=real(Iin);
I2=imag(Iin);
Z2APEX=ZL(1)-LLin;
fprintf(fid_22,'Source  %2.4f *******  *******  %8.2f  %8.2f  %8.2f  %8.2f',Z2APEX,R1,I1,R2,I2);
fprintf(fid_22,'\n');

for ii=1:N,
  R1=real(VA(ii));
  I1=imag(VA(ii));
  R2=real(IA(ii));
  I2=imag(VA(ii));
  fprintf(fid_22,'%6.0f  %2.4f  %2.4f  %2.5f  %8.2f  %8.2f  %8.2f  %8.2f',ii,ZL(ii),L(ii),D(ii),real(VA(ii)),imag(VA(ii)),real(IA(ii)),imag(IA(ii)));
  fprintf(fid_22,'\n');
end

R1=real(Vout);
I1=imag(Vout);
R2=real(Iout);
I2=imag(Iout);
Z2APEX=ZL(N)+LLout;
fprintf(fid_22,'Term.   %2.4f *******  *******  %8.2f  %8.2f  %8.2f  %8.2f',Z2APEX,R1,I1,R2,I2);
fprintf(fid_22,'\nCurrents and voltages relative to  1.0 Watt of input power.');
%fprintf(fid_22,'\n ');
fprintf(fid_22,'\n');
fprintf(fid_22,'\n');
fprintf(fid_22,'Other Design Parameters');
fprintf(fid_22,'\n   Upper Design Frequency (MHz)            :     %6.5f',Fhigh);
fprintf(fid_22,'\n   Lower Design Frequency (MHz)            :     %6.5f',Flow);
fprintf(fid_22,'\n   Tau                                     :     %6.5f',Tau);
fprintf(fid_22,'\n   Sigma                                   :     %6.5f',sigm);
fprintf(fid_22,'\n   Source Impedance (Ohms)                 :     %6.5f',Rs);
fprintf(fid_22,'\n');
%fprintf(fid_22,'\n');
fprintf(fid_22,'   Characteristic Impedance of Source');
fprintf(fid_22,'\n       Transmission Line                   :     %6.5f +j %6.5f',real(ZCin),imag(ZCin));
fprintf(fid_22,'\n   Transminating Impedance (Ohms)          :     %6.5f +j %6.5f',real(Zout),imag(Zout));
fprintf(fid_22,'\n');
%fprintf(fid_22,'\n');
fprintf(fid_22,'   Characteristic Impedance of Termination');
fprintf(fid_22,'\n       Transmission Line                   :     %6.5f +j %6.5f',real(Z0),imag(Z0));
fprintf(fid_22,'\n');
fprintf(fid_22,'\n');
fprintf(fid_22,'Transmission Line of Antenna');
fprintf(fid_22,'\n   Boom Diameter (cm)                      :     %6.5f',DB);
fprintf(fid_22,'\n   Boom Spacing  (cm)                      :     %6.5f',SB);
fprintf(fid_22,'\n   Characteristic Impedance                :     %6.5f +j %6.5f',real(Z0),imag(Z0));
fprintf(fid_22,'\n');
fprintf(fid_22,'\nArray Parameters Analysed at    %6.5f MHz',FMHZ);
fprintf(fid_22,'\n   Array Peak Directivity (dBi)            :     %6.5f',PD0);
fprintf(fid_22,'\n   Array Directivity on Boresight (dBi)    :     %6.5f',D0A);
fprintf(fid_22,'\n   Array Front to Back Ratio (dB)          :     %6.5f',FTBR);

if FSLL>8000
  fprintf(fid_22,'\n   Array-Front-to-Sidelobe Level (dB)      :     *******');
else
  fprintf(fid_22,'\n   Array-Front-to-Sidelobe Level (dB)      :     %6.5f',FSLL);
end
fprintf(fid_22,'\n');
%fprintf(fid_22,'\n');
fprintf(fid_22,'   Input Impedance at shortest element');
fprintf(fid_22,'\n       (Complex - Ohms)                    :     %6.5f +j %6.5f',real(ZinA),imag(ZinA));
fprintf(fid_22,'\n');
fprintf(fid_22,'   Input Impedance at shortest element');
fprintf(fid_22,'\n       (Magnitude - Ohms)                  :     %6.5f ',abs(ZinA));
fprintf(fid_22,'\n   Input Impedance seen by the source  ');
fprintf(fid_22,'\n       (Complex - Ohms)                    :     %6.5f +j %6.5f',real(ZinAS),imag(ZinAS));
fprintf(fid_22,'\n   Input Impedance seen by the source  ');
fprintf(fid_22,'\n       (Magnitude - Ohms)                  :     %6.5f  ',abs(ZinAS));
fprintf(fid_22,'\n   VSWR caused by mismatch between the source');
fprintf(fid_22,'\n       transmission line and the antenna ');
fprintf(fid_22,'\n       and relative to %6.5f +j %6.5f   Ohms :     %6.5f',real(ZCin),imag(ZCin),VSWR);


%***************************************************************
%***************************************************************

function output2(N, ZL, L, D, Z0,design_TITLE, Fhigh, Flow,...
                 Tau, sigm, D0, LD, Rs, LLin,ZCin,...
                  Rin, LLout, Zout, SB, DB)

Alpha = atan((1-Tau)/4.0/sigm)*180/pi;
Z2APEX = LLout + ZL(1);

fid_sum = fopen('LP_des.out','wt');

fprintf(fid_sum,'DIPOLE ARRAY DESIGN');
fprintf(fid_sum,'\n');
fprintf(fid_sum,design_TITLE);
fprintf(fid_sum,'\n');
fprintf(fid_sum,'\n');
fprintf(fid_sum,'Ele.         Z       L       D');
fprintf(fid_sum,'\n');
fprintf(fid_sum,'           (m)     (m)     (cm)');
fprintf(fid_sum,'\n');
  fprintf(fid_sum,'Term.   %2.4f *******  *******',Z2APEX);
fprintf(fid_sum,'\n');
for i=1:N,
  fprintf(fid_sum,'%6.0f  %2.4f  %2.4f  %2.5f',i,ZL(i),L(i),D(i));
  fprintf(fid_sum,'\n');
end
  fprintf(fid_sum,'Source  %2.4f *******  *******',Z2APEX);
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'Design Parameters'); 
fprintf(fid_sum,'\n   Upper Design Frequency (MHz)      :   %6.5f',Fhigh);
fprintf(fid_sum,'\n   Lower Design Frequency (MHz)      :   %6.5f',Flow);
fprintf(fid_sum,'\n   Tau                               :   %6.5f',Tau);
fprintf(fid_sum,'\n   Sigma                             :   %6.5f',sigm);
fprintf(fid_sum,'\n   Alpha (deg)                       :   %6.5f',Alpha);

if D0~=0
  fprintf(fid_sum,'\n   Desired Directivity               :   %6.5f',D0);
else
  fprintf(fid_sum,'\n   Desired Directivity               :   ******');
end

fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'Source and Source Transmission Line'); 
fprintf(fid_sum,'\n   Source Resistance (Ohms)          :   %6.5f',Rs);
fprintf(fid_sum,'\n   Transmission Line Length (m)      :   %6.5f',LLin);
fprintf(fid_sum,'\n   Characteristic Impedance (Ohms)   :   %6.5f +j %6.5f',real(ZCin),imag(ZCin));
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'Antenna and Antenna Transmission Line'); 
fprintf(fid_sum,'\n   Length-to-Diameter Ratio          :   %6.5f',LD);
fprintf(fid_sum,'\n   Boom Diameter (cm)                :   %6.5f',DB);
fprintf(fid_sum,'\n   Boom Spacing  (cm)                :   %6.5f',SB);
fprintf(fid_sum,'\n   Characteristic Impedance (Ohms)   :   %6.5f +j %6.5f',real(Z0),imag(Z0));
fprintf(fid_sum,'\n   Desired Input Impedance (Ohms)    :   %6.5f',Rin);
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'\n'); 
fprintf(fid_sum,'Termination and Termination Transmission Line');
fprintf(fid_sum,'\n   Termination Impedance (Ohms)      :   %6.5f +j %6.5f',real(Zout),imag(Zout));
fprintf(fid_sum,'\n   Transmission Line Length (m)      :   %6.5f',LLout);
fprintf(fid_sum,'\n   Characteristic Impedance (Ohms)   :   %6.5f +j %6.5f',real(Z0),imag(Z0));
 
fclose(fid_sum);

%******************************************************************
%******************************************************************

function [YT]=find_YT(N,DZL,BETA,Z0,Zout)

%******************************************************************
%  This routine calculates the transmission line short circuit
%  admittance matrix, YT.  An approximation is made in certain
%  cases when a division by zero error would normally result.
%******************************************************************

tiny=1.e-12;
big =1.e12;

Y0=1/Z0;

YT=zeros(N,N);

% ELEMENTS 2 TO N-1

if (N>2)

  for ii=2:N-1,
 
    if ( abs( sin( BETA*DZL(ii-1) ) ) <= tiny )
      T1=big;
    else
      T1=1/tan( BETA*DZL(ii-1) );
    end
    
    if ( abs( sin( BETA*DZL(ii) ) ) <= tiny )
      T2=big;
    else
      T2=1/tan(BETA*DZL(ii));
    end

    YT(ii,ii)=-i*(T1+T2)*Y0;

    if ( abs( sin( BETA*DZL(ii-1) ) ) <= tiny )
      T1=sgn( cos( BETA*DZL(ii-1) ) )*big;
    else
      T1=1/sin(BETA*DZL(ii-1));
    end
   
    YT(ii-1,ii)=-i*Y0*T1;

    if ( abs( sin( BETA*DZL(ii) ) ) <= tiny )
      T2=sgn( cos(BETA*DZL(ii) ) )*big;
    else
      T2=1/sin(BETA*DZL(ii));
    end

    YT(ii+1,ii)=-i*Y0*T2;

  end   % end if
 
end     % end if

%ELEMENT N

  if ( abs( sin( BETA*DZL(N-1) ) ) <= tiny )
    T1=big;
  else
    T1=1/tan( BETA*DZL(N-1) );
  end

  if ( abs( Y0*Zout*cos( BETA*DZL(N) )+i*sin( BETA*DZL(N) ) ) == 0 )
    T2=big;
  else
    T2=Y0*(        cos(BETA*DZL(N)) + i*Y0*Zout*sin(BETA*DZL(N)))/...
          (Y0*Zout*cos(BETA*DZL(N)) + i*        sin(BETA*DZL(N)));
  end
    
  if abs(T2)==0
    T2=tiny;
  end
  
  YT(N,N)=-i*Y0*T1+T2;

  if ( abs( sin(BETA*DZL(N-1) ) ) <= tiny )
    T1=sgn(cos(BETA*DZL(N-1)))*big;
  else
    T1=1/sin(BETA*DZL(N-1));
  end

  YT(N-1,N)=-i*Y0*T1;

% ELEMENT 1

  if ( abs( sin( BETA*DZL(1) ) )<= tiny )
    T1=big;
  else
    T1=1/tan(BETA*DZL(1));
  end

  YT(1,1)=-i*Y0*T1;

  if ( abs( sin( BETA*DZL(1) ) ) <= tiny )
    T1=sgn(cos(BETA*DZL(1)))*big;
  else
    T1=1/sin(BETA*DZL(1));
  end

  YT(2,1)=-i*Y0*T1;

%******************************************************************
%******************************************************************

function [ZA]=find_ZA(N,BETA,D,ZL,L)

%******************************************************************
%   This subroutine finds the antenna self and mutual impedances 
%   between elements assuming the current distribution is sinusoidal.
%   Note that 'element' 1 contains the driving current source.
%******************************************************************

ETA    = 120*pi;       % free space imedance
egamma = .5772156649;  % Euler' constant
LAMBDA = 2*pi/BETA;    % wavelength
load trig_int.mat

ZA=zeros(N,N);

for ii=1:N,
  for jj=ii:N,

    LENGTH(1) = L(ii)/2;
    LENGTH(2) = L(jj)/2;
    AELE      = D(ii)/2;
    ZL1       = ZL(ii);
    ZL2       = ZL(jj);

    if (ii==jj)
      BL   = BETA*LENGTH(1)*2;
      BA2L = BETA*AELE^2/LENGTH(1)/2;
      
      x_interp1 = [BL 2*BL];
      y_interp1 = interp1(x_int,Si,x_interp1);
      
      sinint_BL  = y_interp1(1,1); 
      sinint_2BL = y_interp1(1,2); 
      
      x_interp2 = [BL 2*BL 2*BA2L];
      y_interp2 = interp1(x_int,Ci,x_interp2);
      
%      disp('find za 1')
%      keyboard
      
      if (sum(isinf(y_interp2))>0)
         i_inf=find(isinf(y_interp2));
         y_interp2(i_inf)=cosint(x_interp2(i_inf));
      end
      
      if (sum(isnan(y_interp2))>0)
         i_inf=find(isnan(y_interp2));
         y_interp2(i_inf)=cosint(x_interp2(i_inf));
      end
      cosint_BL    = y_interp2(1,1);
      cosint_2BL   = y_interp2(1,2);
      cosint_2BA2L = y_interp2(1,3);
      
 %     disp('find za 2')
 %     keyboard
     
      Rr   = egamma + log(BL) - cosint_BL +...
             .5*sin(BL)*( sinint_2BL - 2*sinint_BL )+...
             .5*cos(BL)*( cosint_2BL - 2*cosint_BL  + egamma + log(BL/2) );

      Xm   =  -cos(BL)*( sinint_2BL - 2*sinint_BL )+...
               sin(BL)*( cosint_2BL - 2*cosint_BL + cosint_2BA2L )+2*sinint_BL;

%      Rr   = egamma + log(BL) - cosint(BL) +...
%             .5*sin(BL)*( sinint(2*BL) - 2*sinint(BL) )+...
%             .5*cos(BL)*( cosint(2*BL) - 2*cosint(BL)  + egamma + log(BL/2) );
%
%      Xm   =  -cos(BL)*( sinint(2*BL) - 2*sinint(BL) )+...
%               sin(BL)*( cosint(2*BL) - 2*cosint(BL) + cosint(2*BA2L) )+2*sinint(BL);

      ZA(ii,jj)= ETA/2/pi*(Rr+i*Xm/2);
      
      ZA(ii,jj)=ZA(ii,jj)/sin(BL/2)^2;

 %   disp('find za 3')
 %    keyboard
      
    else

      DZL   = ZL2 - ZL1;
      L1L2p = LENGTH(1)+LENGTH(2);
      L1L2m = LENGTH(1)-LENGTH(2);
      u0    = BETA*( sqrt(DZL^2 + L1L2p^2) - L1L2p );
      v0    = BETA*( sqrt(DZL^2 + L1L2p^2) + L1L2p );
      u0p   = BETA*( sqrt(DZL^2 + L1L2m^2) - L1L2m );
      v0p   = BETA*( sqrt(DZL^2 + L1L2m^2) + L1L2m );
      u1    = BETA*( sqrt(DZL^2 + LENGTH(1)^2) - LENGTH(1) );
      v1    = BETA*( sqrt(DZL^2 + LENGTH(1)^2) + LENGTH(1) );
      w1    = BETA*( sqrt(DZL^2 + LENGTH(2)^2) + LENGTH(2) );
      y1    = BETA*( sqrt(DZL^2 + LENGTH(2)^2) - LENGTH(2) );
      BLp   = BETA*L1L2p;
      BLm   = BETA*L1L2m;
      BD    = BETA*DZL;

    x_interp1 = [u0 v0];
    x_interp2 = [u1 v1];
    x_interp3 = [w1 y1 BD];
    x_interp4 = [u0p v0p];
    
    y_interp1 = interp1(x_int,Ci,x_interp1);
    y_interp2 = interp1(x_int,Ci,x_interp2);
    y_interp3 = interp1(x_int,Ci,x_interp3);
    y_interp4 = interp1(x_int,Ci,x_interp4);
    
    if sum(isinf(y_interp1))>0
         i_inf=find(isinf(y_interp1));
         y_interp1(i_inf)=cosint(x_interp1(i_inf));
    end
     
    if (sum(isnan(y_interp1))>0)
         i_inf=find(isnan(y_interp1));
         y_interp1(i_inf)=cosint(x_interp1(i_inf));
    end
      
    if sum(isinf(y_interp2))>0
         i_inf=find(isinf(y_interp2));
         y_interp2(i_inf)=cosint(x_interp2(i_inf));
    end
      
    if (sum(isnan(y_interp2))>0)
         i_inf=find(isnan(y_interp2));
         y_interp2(i_inf)=cosint(x_interp2(i_inf));
    end
    
    if sum(isinf(y_interp3))>0
         i_inf=find(isinf(y_interp3));
         y_interp3(i_inf)=cosint(x_interp3(i_inf));
    end
      
    if (sum(isnan(y_interp3))>0)
         i_inf=find(isnan(y_interp3));
         y_interp3(i_inf)=cosint(x_interp3(i_inf));
    end
    
    if sum(isinf(y_interp4))>0
         i_inf=find(isinf(y_interp4));
         y_interp4(i_inf)=cosint(x_interp4(i_inf));
    end
      
    if (sum(isnan(y_interp4))>0)
         i_inf=find(isnan(y_interp4));
         y_interp4(i_inf)=cosint(x_interp4(i_inf));
    end
    y_interp5 = interp1(x_int,Si,x_interp1);
    y_interp6 = interp1(x_int,Si,x_interp2);
    y_interp7 = interp1(x_int,Si,x_interp3);
    y_interp8 = interp1(x_int,Si,x_interp4);
    
    cosint_u0 = y_interp1(1,1);
    cosint_v0 = y_interp1(1,2);
    
    cosint_u1 = y_interp2(1,1);
    cosint_v1 = y_interp2(1,2);
    
    cosint_w1 = y_interp3(1,1);
    cosint_y1 = y_interp3(1,2);
    cosint_BD = y_interp3(1,3);
    
    cosint_u0p = y_interp4(1,1);
    cosint_v0p = y_interp4(1,2);
    
    sinint_u0 = y_interp5(1,1);
    sinint_v0 = y_interp5(1,2);
    
    sinint_u1 = y_interp6(1,1);
    sinint_v1 = y_interp6(1,2);
    
    sinint_w1 = y_interp7(1,1);
    sinint_y1 = y_interp7(1,2);
    sinint_BD = y_interp7(1,3);
    
    sinint_u0p = y_interp8(1,1);
    sinint_v0p = y_interp8(1,2);
    
      Rr    = cos(BLp)*( cosint_u0 + cosint_v0 - cosint_u1 - cosint_v1-...
                         cosint_w1 - cosint_y1 + 2*cosint_BD)+...
              cos(BLm)*( cosint_u0p+ cosint_v0p- cosint_u1 - cosint_v1-...
                         cosint_w1 - cosint_y1 + 2*cosint_BD)+...
              sin(BLp)*(-sinint_u0 + sinint_v0+ sinint_u1 - sinint_v1-...
                         sinint_w1 + sinint_y1)+...
              sin(BLm)*(-sinint_u0p+ sinint_v0p+ sinint_u1 - sinint_v1+...
                         sinint_w1 - sinint_y1);
     
      Xm   = cos(BLp)*(-sinint_u0 - sinint_v0 + sinint_u1 + sinint_v1 +...
                        sinint_w1 + sinint_y1 - 2*sinint_BD)+...
             cos(BLm)*(-sinint_u0p - sinint_v0p + sinint_u1 + sinint_v1+...
                        sinint_w1 + sinint_y1 - 2*sinint_BD)+...
             sin(BLp)*(-cosint_u0 + cosint_v0 + cosint_u1 - cosint_v1 -...
                        cosint_w1 + cosint_y1 )+...
             sin(BLm)*(-cosint_u0p + cosint_v0p + cosint_u1 - cosint_v1 +...
                        cosint_w1 - cosint_y1 );



%      Rr    = cos(BLp)*( cosint(u0) + cosint(v0) - cosint(u1) - cosint(v1)-...
%                         cosint(w1) - cosint(y1) + 2*cosint(BD))+...
%              cos(BLm)*( cosint(u0p)+ cosint(v0p)- cosint(u1) - cosint(v1)-...
%                         cosint(w1) - cosint(y1) + 2*cosint(BD))+...
%              sin(BLp)*(-sinint(u0) + sinint(v0)+ sinint(u1) - sinint(v1)-...
%                         sinint(w1) + sinint(y1))+...
%              sin(BLm)*(-sinint(u0p)+ sinint(v0p)+ sinint(u1) - sinint(v1)+...
%                         sinint(w1) - sinint(y1));
%     
%      Xm   = cos(BLp)*(-sinint(u0) - sinint(v0) + sinint(u1) +sinint(v1)+...
%                        sinint(w1) + sinint(y1) - 2*sinint(BD))+...
%             cos(BLm)*(-sinint(u0p)-sinint(v0p) + sinint(u1) + sinint(v1)+...
%                        sinint(w1) + sinint(y1) - 2*sinint(BD))+...
%             sin(BLp)*(-cosint(u0)+cosint(v0)+cosint(u1)-cosint(v1)-...
%                        cosint(w1)+cosint(y1))+...
%             sin(BLm)*(-cosint(u0p)+cosint(v0p)+cosint(u1)-cosint(v1)+...
%                        cosint(w1)-cosint(y1));

      ZA(ii,jj)=ETA/4/pi*(Rr+i*Xm);

      ZA(ii,jj)=ZA(ii,jj)/sin(BETA*LENGTH(1))/sin(BETA*LENGTH(2));

      ZA(jj,ii)=ZA(ii,jj);

    end
  end
end

%******************************************************************
%******************************************************************

function GAIN=pattern(BETA,Phi,Theta,L,ZL,N,Iel,RADIOUS)

ms=30;
TINY=1.e-12;

ETA=120*pi;

KB(1)=sin( Theta*pi/180 )*cos( Phi*pi/180 );
KB(2)=sin( Theta*pi/180 )*sin( Phi*pi/180 );
KB(3)=cos( Theta*pi/180 );

PAT=0;

for ii=1:N,
  PHASE  = mod( BETA*( RADIOUS-ZL(ii)*KB(3) ),2*pi );
  IsoPAT = i*ETA*exp(-i*PHASE)/2/pi/RADIOUS;

  KL2   = BETA * L(ii)/2;
  cosTH1= KB(2);
  sinTH1= sqrt(1-cosTH1^2);

  if (abs(sinTH1)>TINY)
    DipPat = IsoPAT*(cos(KL2*cosTH1) - cos(KL2))/sinTH1;
  else
    DipPat=0;
  end

  PAT=PAT+Iel(ii)*DipPat;
end
  
  if abs(PAT)>0
    GAIN=10*log10(abs(PAT)^2/2/ETA);
  else
    GAIN=-1000;
  end








