%***********************************************************************
%
%     SECTORAL HORN: UNBOUNDED MEDIUM
%
%     PROGRAM AUTHOR--  WILLIAM V. ANDREW
%                       DEPARTMENT OF ELECTRICAL ENGINEERING
%                       ARIZONA STATE UNIVERSITY
%                       TEMPE, ARIZONA  85287-7206
%                       (602) 965-5311
%
%     DATE OF THIS VERSION--  Nov. 3, 1995
%
%     This MATLAB M-file will produce the FDTD solution
%     of a sectoral (2-D) Perfectly Electric Conducting
%     (PEC) horn antenna excited by a sinusoidal voltage
%     in a TEz computational domain. The computational
%     domain is truncated with a Berenger Perfectly Matched
%     Layer (PML) absorbing boundary condition whose depth
%     in layers is set by the variable NPMLS. The PML is
%     introduced to eliminate reflections from the grid
%     truncation and to simulate an outgoing traveling
%     wave propagating in an unbounded medium. The M-file
%     can also create a movie.  For example, you can create
%     a movie which is 70 frames long by taking a picture of 
%     the computational domain every 3rd time step.
%     
%     To execute this M-file, type ``te_horn'' at the
%     MATLAB prompt. The file will save each frame of
%     the movie and then write the entire movie to a file
%     named ``te_horn.mat''.
%
%     To play the movie at any time after it has been 
%     created and saved to file ``te_horn.mat'' just
%     execute the following MATLAB commands:
%
%          load te_horn.mat
%          movie(M,n,fps) 
%
%     where n is the number of times to play the movie
%     and fps is the number of frames per second. 
%
%     This M-file will not work with the Student edition
%     of MATLAB due to the restrictions on array size.
%     Therefore, this M-file will work only with the
%     Professional edition of MATLAB. The movie which
%     this file creates is approximately 10.6 Mbytes in
%     size. Therefore, the available RAM memory or the
%     swap space on whatever operating system must be
%     large enough to accommodate a file this large.
%
%     The horn is modeled by setting the necessary
%     FDTD update equation coefficients to represent
%     PEC material (sigma=infinity). 
%     The cell size of the space is:
%         dx = 0.0025 meters
%     The time step is:
%         dt = 4.23e-12 seconds
%     The frequency of excitation is:
%         freq = 9.84252 GHz
%     The wavelength is:
%         lambda = 12*dx = 0.0305 meters
%
%     The flare section of the horn is staircased. As
%     modeled, the horn looks like:
%
%             jc-7                        jc+7                  
%               |                           |
%             \ | /                       \ | / 
%              \|/                         \|/
%               `                           `   
%
%   ic+11       |_                         _| 
%   ic+10         |                       |
%   ic+9          |_                     _|
%   ic+8            |                   |
%   ic+7            |_                 _|
%   ic+6              |               |
%   ic+5              |_             _|
%   ic+4                |           |
%   ic+3 ------------>  |_         _|
%   ic+2                  |       |
%   ic+1                  |       |
%   ic                    |       |
%   ic-1                  |       |
%   ic-2                  |       |
%   ic-3                  |       |
%   ic-4                  |       |     `|'  <--- Ex Field Component
%   ic-5                  |       |
%   ic-6                  |       |     `-'  <--- Ey Field Component
%   ic-7                  |       |
%   ic-8                  |       |
%   ic-9                  |       |   
%   ic-10-------------->  |_ _ _ _|  <----- Excitation Plane
%   ic-11                 |       |  
%   ic-12                 |       |
%   ic-13---------------> |_ _ _ _|
%
%                         .       .
%                        /|\     /|\
%                       / | \   / | \
%                         |       |
%                       jc-2    jc+2
%
%***********************************************************************

clear      

REPEAT = input('How many times to repeat? ');
FPS    = input('Enter frames per second... ');

%***********************************************************************
%     Initialize some constants
%***********************************************************************

      npmls=8;                 % Depth of PML region in # of cells
      
      nmax=210;                % Number of time steps

      ie=100;                 
      ib=ie+1;
      ic=ie/2-20;
      ip=ie-npmls;
      
      je=100;
      jb=je+1;
      jc=je/2;
      jp=je-npmls;

      pi=4.0*atan(1.0);
      muo=4.0*pi*1.0e-7;       % Permeability of free space
      epso=8.854e-12;          % Permittivity of free space
      co=1.0/sqrt(muo*epso);   % Speed of light in free space
      aimp=sqrt(muo/epso);     % Wave impedance in free space
      freq=9.84252e+09;        % Frequency of excitation
      lambda=co/freq;          % Wavelength of excitation
      dx=lambda/12.0;          % FDTD cell size
      dt=dx/co/2.0;            % Time step size

%***********************************************************************
%     .... Set up the Berenger PML ABC material constants ....
%***********************************************************************

      sigmax=-3.0*epso*co*log(1e-5)/(2.0*dx*npmls);
      rhomax=sigmax*(aimp^2);

      for m=1:npmls;
       sig(m)=sigmax*((m-0.5)/(npmls+0.5))^2;
       rho(m)=rhomax*(m/(npmls+0.5))^2;
      end;

%***********************************************************************
%     .... Set up constants needed in the FDTD equations for the ....
%     .... Berenger PML ABCs (exponential difference expressions)....
%***********************************************************************

      for m=1:npmls;
       re=sig(m)*dt/epso;
       rm=rho(m)*dt/muo;
       ca(m)=exp(-re);
       cb(m)=-(exp(-re)-1.0)/sig(m)/dx;
       da(m)=exp(-rm);
       db(m)=-(exp(-rm)-1.0)/rho(m)/dx;
      end;

%***********************************************************************
%   Initialize all of the matrices for the field components HZ, HZX,
%   HZY, EX, EY, CAEX, CAEY, DAHZX, DAHZY, CBEX, CBEY, DBHZX, and DBHZY.
%***********************************************************************

      for i=1:ie;
       for j=1:jb;
	ex(i,j)=0.0;
	caex(i,j)=1.0;           % Free space
	cbex(i,j)=dt/epso/dx;    % Free space
       end;
      end;

      for i=1:ib;
       for j=1:je;
	ey(i,j)=0.0;
	caey(i,j)=1.0;          % Free space
	cbey(i,j)=dt/epso/dx;   % Free space
       end;
      end;

      for i=1:ie;
       for j=1:je;
	hz(i,j)=0.0;
	hzx(i,j)=0.0;
	dahzx(i,j)=1.0;           % Free space
	dbhzx(i,j)=dt/muo/dx;     % Free space
	hzy(i,j)=0.0;
	dahzy(i,j)=1.0;           % Free space
	dbhzy(i,j)=dt/muo/dx;     % Free space
       end;
      end;

%*******************************************************************
%     Initialize all of the matrices for the Berenger PML absorbing
%     boundaries.
%*******************************************************************

%<<<<<<<<<<<<<<<<<<<<<<<<< Ex Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Left and Right PML Regions .....

      for i=2:ie;
       for j=2:npmls+1;
	m=npmls+2-j;
	caex(i,j)=ca(m);
	cbex(i,j)=cb(m);
       end;
       for j=jp+1:je;
	m=j-jp;
	caex(i,j)=ca(m);
	cbex(i,j)=cb(m);
       end;
      end;

%<<<<<<<<<<<<<<<<<<<<<<<<< Ey Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Back and Front PML Regions .....

      for j=2:je;
       for i=2:npmls+1;
	m=npmls+2-i;
	caey(i,j)=ca(m);
	cbey(i,j)=cb(m);
       end;
       for i=ip+1:ie;
	m=i-ip;
	caey(i,j)=ca(m);
	cbey(i,j)=cb(m);
       end;
      end;

%<<<<<<<<<<<<<<<<<<<<<<<<< Hz Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Left and Right PML Regions .....

      for i=2:ie;
       for j=1:npmls;
	m=npmls+1-j;
	dahzy(i,j)=da(m);
	dbhzy(i,j)=db(m);
       end;
       for j=jp+1:je;
	m=j-jp;
	dahzy(i,j)=da(m);
	dbhzy(i,j)=db(m);
       end;
      end;

%        ..... Front and Back PML Regions .....

      for j=2:je;
       for i=1:npmls;
	m=npmls+1-i;
	dahzx(i,j)=da(m);
	dbhzx(i,j)=db(m);
       end;
       for i=ip+1:ie;
	m=i-ip;
	dahzx(i,j)=da(m);
	dbhzx(i,j)=db(m);
       end;
      end;

%***********************************************************************
%          ..... Material coefficients PEC of horn antenna .....
%***********************************************************************

      for i=ic-13:ic+2;
        caex(i,jc-2)=-1.0;
        cbex(i,jc-2)= 0.0;
        caex(i,jc+2)=-1.0;
        cbex(i,jc+2)= 0.0;
      end;

      caex(ic+3,jc-3)=-1.0;
      cbex(ic+3,jc-3)= 0.0;

      caex(ic+3,jc+3)=-1.0;
      cbex(ic+3,jc+3)= 0.0;

      caex(ic+4,jc-3)=-1.0;
      cbex(ic+4,jc-3)= 0.0;

      caex(ic+4,jc+3)=-1.0;
      cbex(ic+4,jc+3)= 0.0;

      caex(ic+5,jc-4)=-1.0;
      cbex(ic+5,jc-4)= 0.0;

      caex(ic+5,jc+4)=-1.0;
      cbex(ic+5,jc+4)= 0.0;

      caex(ic+6,jc-4)=-1.0;
      cbex(ic+6,jc-4)= 0.0;

      caex(ic+6,jc+4)=-1.0;
      cbex(ic+6,jc+4)= 0.0;

      caex(ic+7,jc-5)=-1.0;
      cbex(ic+7,jc-5)= 0.0;

      caex(ic+7,jc+5)=-1.0;
      cbex(ic+7,jc+5)= 0.0;

      caex(ic+8,jc-5)=-1.0;
      cbex(ic+8,jc-5)= 0.0;

      caex(ic+8,jc+5)=-1.0;
      cbex(ic+8,jc+5)= 0.0;

      caex(ic+9,jc-6)=-1.0;
      cbex(ic+9,jc-6)= 0.0;

      caex(ic+9,jc+6)=-1.0;
      cbex(ic+9,jc+6)= 0.0;

      caex(ic+10,jc-6)=-1.0;
      cbex(ic+10,jc-6)= 0.0;

      caex(ic+10,jc+6)=-1.0;
      cbex(ic+10,jc+6)= 0.0;

      caex(ic+11,jc-7)=-1.0;
      cbex(ic+11,jc-7)= 0.0;

      caex(ic+11,jc+7)=-1.0;
      cbex(ic+11,jc+7)= 0.0;

      for j=jc-2:jc+1;
       caey(ic-13,j)=-1.0;
       cbey(ic-13,j)= 0.0;
      end;

      caey(ic+3,jc-3)=-1.0;
      cbey(ic+3,jc-3)= 0.0;

      caey(ic+3,jc+2)=-1.0;
      cbey(ic+3,jc+2)= 0.0;

      caey(ic+5,jc-4)=-1.0;
      cbey(ic+5,jc-4)= 0.0;

      caey(ic+5,jc+3)=-1.0;
      cbey(ic+5,jc+3)= 0.0;

      caey(ic+7,jc-5)=-1.0;
      cbey(ic+7,jc-5)= 0.0;

      caey(ic+7,jc+4)=-1.0;
      cbey(ic+7,jc+4)= 0.0;

      caey(ic+9,jc-6)=-1.0;
      cbey(ic+9,jc-6)= 0.0;

      caey(ic+9,jc+5)=-1.0;
      cbey(ic+9,jc+5)= 0.0;

      caey(ic+11,jc-7)=-1.0;
      cbey(ic+11,jc-7)= 0.0;

      caey(ic+11,jc+6)=-1.0;
      cbey(ic+11,jc+6)= 0.0;

%***********************************************************************
%          .....TIME-STEPPING LOOP.....
%***********************************************************************

      for n=1:nmax;

%***********************************************************************
%          .....EX FIELD UPDATE.....
%***********************************************************************

       ex(1:ie,2:je)=caex(1:ie,2:je).*ex(1:ie,2:je)+...
	cbex(1:ie,2:je).*(hz(1:ie,2:je)-hz(1:ie,1:je-1));

%***********************************************************************
%          .....EY FIELD UPDATE.....
%***********************************************************************

      ey(2:ie,1:je)=caey(2:ie,1:je).*ey(2:ie,1:je)+...
       cbey(2:ie,1:je).*(hz(1:ie-1,1:je)-hz(2:ie,1:je));
      
%***********************************************************************
%     ..... Hard Source ramped sinusoidal excitation .....
%***********************************************************************

      for j=jc-2:jc+1;
        ey(ic-10,j)=(1.0-exp(-((n/20.0)^2)))*...
                    aimp*sin(2.0*pi*freq*n*dt);
      end;

%***********************************************************************
%          .....HZ FIELD UPDATE.....            
%***********************************************************************

      hzx(1:ie,1:je)=dahzx(1:ie,1:je).*hzx(1:ie,1:je)+...
       dbhzx(1:ie,1:je).*(ey(1:ie,1:je)-ey(2:ib,1:je));

      hzy(1:ie,1:je)=dahzy(1:ie,1:je).*hzy(1:ie,1:je)+...
       dbhzy(1:ie,1:je).*(ex(1:ie,2:jb)-ex(1:ie,1:je));

      hz(1:ie,1:je)=hzx(1:ie,1:je)+hzy(1:ie,1:je);

%***********************************************************************
%    .....Create the movie frame by frame.....
%    .....Take a frame every 3rd time step.....
%***********************************************************************

      if rem(n,3)==0;
       s=int2str(n);
       n2=n/3;
       clf;
       pcolor(log10(abs(ey+0.000001)));
       axis([0 100 0 100]);
       caxis([-6 3]);
       shading interp;
       if n==3;
	M=moviein(70);
       end;
       t2=['TEz 2D Horn Antenna. Time step #',s];
       title(t2);
       hold;
       M(:,n2)=getframe;
      end;

%***********************************************************************
%     End time step loop
%***********************************************************************

      end;

%***********************************************************************
%     Save the movie to file ``te_horn.mat''
%***********************************************************************

      save te_horn M;

%***********************************************************************
%     Replay the movie 5 times at 7 frames per second
%***********************************************************************

      movie(M,REPEAT,FPS)







