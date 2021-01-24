%***********************************************************************
%
%     LINE SOURCE EXCITED WITH A TIME DERIVATIVE GAUSSIAN PULSE:
%     UNBOUNDED MEDIUM
%
%     PROGRAM AUTHOR--  WILLIAM V. ANDREW
%                       DEPARTMENT OF ELECTRICAL ENGINEERING
%                       ARIZONA STATE UNIVERSITY
%                       TEMPE, ARIZONA  85287-7206
%                       (602) 965-5311
%
%     DATE OF THIS VERSION--  Nov. 3, 1995
%
%     This MATLAB M-file will produce the FDTD solution of a z-directed
%     line source in a two-dimensional TMz computational domain excited
%     by a time derivative Gaussian pulse. The computational domain
%     is truncated with a Berenger Perfectly Matched Layer (PML)
%     absorbing boundary condition whose depth in layers is set by the
%     variable NPMLS. The PML is introduced to eliminate reflections
%     from the grid truncation and to simulate an outgoing traveling
%     wave propagating in an unbounded medium. The M-file
%     can also create a movie.  For example, you can create
%     a movie which is 70 frames long by taking a picture of 
%     the computational domain every 3rd time step.
%     
%     To execute the M-file, type ``tm_open'' at the MATLAB prompt.
%     The file will save each frame of the movie, replay the movie 5
%     times at 7 frames per second and then write the entire movie
%     to a file named ``tm_open.mat''.
%
%     To play the movie at any time after it has been created and
%     saved to the file ``tm_open.mat'' just execute the following
%     MATLAB commands:
%
%          load tm_open.mat
%          movie(M,n,fps) 
%
%     where n is the number of times to play the movie and 
%     fps is the number of frames per second. 
%
%     This M-file will not work with the Student edition of MATLAB
%     due to the restrictions on array size. Therefore, this M-file
%     will work only with the Professional edition of MATLAB. The
%     movie which this file creates is approximately 5.6 Mbytes in size.
%     Therefore, the available RAM memory or the swap space on whatever
%     operating system must be large enough to accommodate a file 
%     this large.
%
%***********************************************************************

clear      

REPEAT = input('How many times to repeat? ');
FPS    = input('Enter frames per second... ');

%***********************************************************************
%     Initialize some constants
%***********************************************************************

      npmls=6;                 % Depth of PML region in # of cells
      
      nmax=111;                % Number of time steps

      ie=50;                 
      ib=ie+1;
      ic=ie/2+1;
      ip=ie-npmls;
      
      je=50;
      jb=je+1;
      jc=je/2+1;
      jp=je-npmls;

      pi=4.0*atan(1.0);
      muo=4.0*pi*1.0e-7;       % Permeability of free space
      epso=8.854e-12;          % Permittivity of free space
      co=1.0/sqrt(muo*epso);   % Speed of light in free space
      aimp=sqrt(muo/epso);     % Wave impedance in free space
%      freq=9.84252e+09;
%      lambda=co/freq;
      dx=0.003;                % FDTD cell size
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
%     Initialize all of the matrices for the field components HX, HY,
%     EZX, EZY, DAHX, DAHY, CAEZX, CAEZY, DBHX, DBHY, CBEZX, and CBEZY.
%***********************************************************************

      for i=1:ib;
       for j=1:jb;
	ez(i,j)=0.0;
	ezx(i,j)=0.0;
	caezx(i,j)=1.0;          % Free space
	cbezx(i,j)=dt/epso/dx;   % Free space
	ezy(i,j)=0.0;
	caezy(i,j)=1.0;          % Free space
	cbezy(i,j)=dt/epso/dx;   % Free space
       end;
      end;

      for i=1:ib;
       for j=1:je;
	hx(i,j)=0.0;
	dahx(i,j)=1.0;           % Free space
	dbhx(i,j)=dt/muo/dx;     % Free space
       end;
      end;

      for i=1:ie;
       for j=1:jb;
	hy(i,j)=0.0;
	dahy(i,j)=1.0;           % Free space
	dbhy(i,j)=dt/muo/dx;     % Free space
       end;
      end;

%*******************************************************************
%     Initialize all of the matrices for the Berenger PML absorbing
%     boundaries.
%*******************************************************************

%<<<<<<<<<<<<<<<<<<<<<<<<< Ez Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Left and Right PML Regions .....

      for i=2:ie;
       for j=2:npmls+1;
	m=npmls+2-j;
	caezy(i,j)=ca(m);
	cbezy(i,j)=cb(m);
       end;
       for j=jp+1:je;
	m=j-jp;
	caezy(i,j)=ca(m);
	cbezy(i,j)=cb(m);
       end;
      end;

%        ..... Back and Front PML Regions .....

      for j=2:je;
       for i=2:npmls+1;
	m=npmls+2-i;
	caezx(i,j)=ca(m);
	cbezx(i,j)=cb(m);
       end;
       for i=ip+1:ie;
	m=i-ip;
	caezx(i,j)=ca(m);
	cbezx(i,j)=cb(m);
       end;
      end;

%<<<<<<<<<<<<<<<<<<<<<<<<< Hx Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Left and Right PML Regions .....

      for i=2:ie;
       for j=1:npmls;
	m=npmls+1-j;
	dahx(i,j)=da(m);
	dbhx(i,j)=db(m);
       end;
       for j=jp+1:je;
	m=j-jp;
	dahx(i,j)=da(m);
	dbhx(i,j)=db(m);
       end;
      end;

%<<<<<<<<<<<<<<<<<<<<<<<<< Hy Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%        ..... Front and Back PML Regions .....

      for j=2:je;
       for i=1:npmls;
	m=npmls+1-i;
	dahy(i,j)=da(m);
	dbhy(i,j)=db(m);
       end;
       for i=ip+1:ie;
	m=i-ip;
	dahy(i,j)=da(m);
	dbhy(i,j)=db(m);
       end;
      end;

%***********************************************************************
%          .....TIME-STEPPING LOOP.....
%***********************************************************************

      for n=1:nmax;

%***********************************************************************
%  .....Set up time derivative Gaussian pulse excitation voltage.....
%***********************************************************************

      b = 25.0;
      dum = 4.0/b/dt*(n*dt-b*dt);
      voltage = 2.0*dum*exp(-(dum^2));

%***********************************************************************
%          .....EZ FIELD UPDATE.....
%***********************************************************************

       ezx(2:ie,2:je)=caezx(2:ie,2:je).*ezx(2:ie,2:je)+...
	cbezx(2:ie,2:je).*(hy(2:ie,2:je)-hy(1:ie-1,2:je));

       ezy(2:ie,2:je)=caezy(2:ie,2:je).*ezy(2:ie,2:je)+...
	cbezy(2:ie,2:je).*(hx(2:ie,1:je-1)-hx(2:ie,2:je));

       ez(2:ie,2:je)=ezx(2:ie,2:je)+ezy(2:ie,2:je);

%***********************************************************************
%     ..... Hard Source excitation .....
%***********************************************************************

       ez(ic,jc)=voltage/dx;

%***********************************************************************
%          .....HX FIELD UPDATE.....            
%***********************************************************************

      hx(1:ib,1:je)=dahx(1:ib,1:je).*hx(1:ib,1:je)+...
       dbhx(1:ib,1:je).*(ez(1:ib,1:je)-ez(1:ib,2:jb));

%***********************************************************************
%          .....HY FIELD UPDATE.....
%***********************************************************************

      hy(1:ie,1:jb)=dahy(1:ie,1:jb).*hy(1:ie,1:jb)+...
       dbhy(1:ie,1:jb).*(ez(2:ib,1:jb)-ez(1:ie,1:jb));

%***********************************************************************
%    .....Create the movie frame by frame.....
%    .....Take a frame every 3rd time step.....
%***********************************************************************

      if rem(n,3)==0;
       s=int2str(n);
       n2=n/3;
       clf;
       pcolor(ez);
       axis([0 50 0 50]);
       caxis([-50 50]);
       shading interp;
       if n==3;
	M=moviein(37);
       end;
       t2=['TMz Gaussian Derivative Pulse. Time step #',s];
       title(t2);
       hold;
       M(:,n2)=getframe;
      end;

%***********************************************************************
%     End time step loop
%***********************************************************************

      end;

%***********************************************************************
%     Save the movie to file ``tm_open.mat''
%***********************************************************************

      save tm_open M;

%***********************************************************************
%     Replay the movie 5 times at 7 frames per second
%***********************************************************************

      movie(M,REPEAT,FPS)







