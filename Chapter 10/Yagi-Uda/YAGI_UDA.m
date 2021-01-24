%***********************************************************************
%   YAGI_UDA.M    
%******************************************************************
%     THIS IS A MATLAB M-FILE THAT COMPUTES, FOR THE YAGI-UDA ARRAY, 
%     THE 
%     
%       I.    FAR-ZONE E- AND H-PLANE AMPLITUDE PATTERNS (IN dB)
%       II.   DIRECTIVITY OF THE ARRAY (IN dB)
%       III.  E-PLANE HALF-POWER BEAMWIDTH (IN DEGREES)
%       IV.   H-PLANE HALF-POWER BEAMWIDTH (IN DEGREES)
%       V.    E-PLANE FRONT-TO-BACK RATIO (IN dB)
%       VI.   H-PLANE FRONT-TO-BACK RATIO (IN dB)
%
%     THE PROGRAM IS BASED ON POCKLINGTON'S INTEGRAL EQUATION 
%     FORMULATION OF SECTION 10.3.3, EQUATIONS (10-42) - (10-65a).  
%     M ENTIRE DOMAIN COSINUSOIDAL (FOURIER) BASIS MODES ARE USED 
%     FOR EACH OF THE ANTENNA ELEMENTS.
%
%       ** INPUT PARAMETERS BY USER
%       1. M     =   NUMBER OF ENTIRE DOMAIN BASIS MODES
%       2. N     =   NUMBER OF ANTENNA ELEMENTS
%       3. L     =   LENGTH OF EACH ELEMENT (IN WAVELENGTHS)
%       4. ALPHA =   RADIUS OF EACH ELEMENT (IN WAVELENGTHS)
%       5. S     =   SEPARATION BERWEEN THE ELEMENTS (IN WAVELENGTHS)
%
%       ** NOTES
%       1.  REFER TO FIGURE 10.19 FOR THE GEOMETRY. 
%       2.  DRIVER ELEMENT IS LOCATED AT THE ORIGIN.
%       3.  FIRST ELEMENT (N=1) IS THE FIRST DIRECTOR.
%       4.  REFLECTOR IS THE N-1 ELEMENT; ONLY ONE REFLECTOR. 
%       5.  DRIVEN ELEMENT IS N.
%
%     THE FORMULATION OF THE PROBLEM IS BASED ON THE PAPER `ANALYSIS OF
%     YAGI-UDA-TYPE ANTENNAS' BY GARY A. THIELE, IEEE TRANS. ANTENNAS
%     PROPAGAT., VOL. 17, JAN. 1969.  
%     ******************************************************************
%     Written by: Mingwei Hsu, Arizona State University
%
%     ******************************************************************

function [] = yagi_uda

%     Declare global variables
global MMAX NMAX Z RHO N2 NMODE L
MMAX = 30;
NMAX = 30;

close all;

%
%     Choice of output
%
fprintf (1, '\n\n');
fprintf (1, '   OUTPUT DEVICE OPTION\n');
fprintf (1, '      OPTION (1): SCREEN\n');
fprintf (1, '      OPTION (2): OUTPUT FILE\n\n');
device = input ('   OUTPUT DEVICE = ', 's');
device = str2num (device);

if (device == 1)
    fid = 1;
elseif (device == 2)
    filename = input ('\n   INPUT THE DESIRED OUTPUT FILENAME = ', 's');
    fid = fopen (filename, 'wt');
else
    fprintf (1, '\n');
    fprintf (1, '   ***ERROR***\n');
    fprintf (1, '   OUTPUT DEVICE NUMBER SHOULD BE EITHER 1 OR 2\n\n\n');
    return
end

%
% INPUT THE LENGTH OF THE ELEMENTS, L, AND THEIR RELATIVE SEPARATION IN FREE
% SPACE WAVELENGTHS, S.  THE VARIABLE YP DEFINES THE ABSOLUTE DISTANCE OF
% THE ELEMENTS ALONG THE Y-AXIS, WITH THE DRIVEN ELEMENT AT THE ORIGIN.      
% ALPHA IS THE ELEMENT WIRE RADIUS IN WAVELENTGHS.
%
%     ***********************************************************************
%     INPUT NUMBER OF MODES PER ELEMENT
%     ---------------------------------
M = input ('\n   NUMBER OF MODES PER ELEMENT (A POSITIVE INTEGER) = ', 's');

M = str2num (M);
M = round (M);
if (M > MMAX)
    fprintf (1, '\n   *** ERROR: Need to increase MMAX in the program.\n\n');
    return
elseif (M <= 0)
    fprintf (1, '\n   *** ERROR: The number has to be greater than 0!\n\n');
    return
elseif isempty (M)   % If the user enters a value other than a number
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return   
end

%     INPUT NUMBER OF ELEMENTS
%     ------------------------
N = input ('   NUMBER OF ELEMENTS (A POSITIVE INTEGER GREATER THAN 1) = ', 's');

N = str2num (N);
N = round (N);
if (N > NMAX)
    fprintf (1, '\n   *** ERROR: Need to increase NMAX in the program.\n\n');
    return
elseif (N <= 0)
    fprintf (1, '\n   *** ERROR: The number has to be greater than 0!\n\n');
    return
elseif isempty (N)   % If the user enters a value other than a number
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end

fprintf (1, '\n');

%     INPUT ELEMENT LENGTHS AND ELEMENT SEPARATION DISTANCES
%     ------------------------------------------------------
%     ---> ELEMENT LENGTHS (IN WAVELENGTHS)
if (N > 3)
    fprintf (1, '   DO ALL DIRECTORS HAVE THE SAME LENGTH?\n');
    ANS = input ('   ANSWER: (Y OR N) ...... ', 's');
else
    ANS = 'N';
end

fprintf (1, '\n');

%     THE LENGTH OF THE DIRECTORS
if (ANS == 'Y') | (ANS == 'y')
    LDIR = input ('   THE UNIFORM LENGTH (in WAVELENGTHS) OF THE DIRECTOR = ', 's');
    LDIR = str2num (LDIR);
    if (isempty (LDIR)) | (LDIR <= 0)
        fprintf (1, '\n   ***** ERROR *****\n');
        fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
        return      
    end
    L = LDIR * ones (1, N-2);
    
elseif (ANS == 'N') | (ANS == 'n')
    a = 1;
    while a <= (N-2)
        fprintf (1, '   LENGTH (in WAVELENGTHS) OF DIRECTOR # %2d =', a);
        b = input (' ', 's');
        b = str2num (b);
        if (isempty (b)) | (b <= 0)
            fprintf (1, '\n   ***** ERROR *****\n');
            fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
            return
        else
            L (a) = b;
        end
        a = a + 1;
    end
else
    
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return
end

%     GET THE LENGTH OF THE REFLECTOR
b = input ('   LENGTH (in WAVELENGTHS) OF THE REFLECTOR = ', 's');
b = str2num (b);
if (isempty (b)) | (b <= 0)
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end
L (N-1) = b;

%     GET THE LENGTH OF THE DRIVEN ELEMENT
b = input ('   LENGTH (in WAVELENGTHS) OF THE DRIVEN ELEMENT = ', 's');
b = str2num (b);
if (isempty (b)) | (b <= 0)
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end
L (N) = b;

%     INPUT ELEMENT SEPARATION DISTANCES
%     ------------------------------------------------------

%     ELEMENT SEPARATION BETWEEN THE DRIVEN ELEMENT AND THE 1ST DIRECTOR
b = input ('\n   SEPARATION (in WAVELENGTHS) BETWEEN DRIVEN ELEMENT & 1ST DIRECTOR = ', 's');
b = str2num (b);
if (isempty (b)) | (b <= 0)
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end

S_1 = b;

if (N > 3)
    fprintf (1, '\n   IS THE SEPARATION BETWEEN DIRECTORS UNIFORM?\n');
    ANS = input ('   ANSWER: (Y OR N) ...... ', 's');
    fprintf (1, '\n');
else
    ANS = 'N';
end

%     THE SEPARATION DISTANCES OF THE DIRECTORS
if (ANS == 'Y') | (ANS == 'y')
    SDIR = input ('   THE UNIFORM SEPARATION (in WAVELENGTHS) BETWEEN DIRECTORS = ', 's');
    SDIR = str2num (SDIR);
    if (isempty (SDIR)) | (SDIR <= 0)
        fprintf (1, '\n   ***** ERROR *****\n');
        fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
        return      
    end
    S = SDIR * ones (1, N-2);
    
elseif (ANS == 'N') | (ANS == 'n')
    a = 2;
    while a <= (N-2)
        fprintf (1, '   SEPARATION (in WAVELENGTHS) BETWEEN DIRECTORS # %2d AND # %2d =', a-1, a);
        b = input (' ', 's');
        b = str2num (b);
        if (isempty (b)) | (b <= 0)
            fprintf (1, '\n   ***** ERROR *****\n');
            fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
            return
        else
            S (a) = b;
        end
        a = a + 1;
    end
    
else
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return
end

S (1) = S_1;

%     ELEMENT SEPARATION BETWEEN THE DRIVEN ELEMENT AND THE REFLECTOR
b = input ('   SEPARATION (in WAVELENGTHS) BETWEEN REFLECTOR & DRIVEN ELEMENT = ', 's');
b = str2num (b);
if (isempty (b)) | (b <= 0)
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end
S (N-1) = b;

%     RADIUS OF EACH ELEMENT
b = input ('\n   RADIUS (in WAVELENGTHS) FOR ALL ELEMENTS USED = ', 's');
b = str2num (b);
if (isempty (b)) | (b <= 0)
    fprintf (1, '\n   ***** ERROR *****\n');
    fprintf (1, '   INPUT DATA ARE NOT OF THE RIGHT FORMAT.\n\n');
    return      
end
ALPHA = b;


%     ************************************************************************************
%     ************************************************************************************

%     Echo all the input parameters
%     -----------------------------
fprintf (fid, '\n   ******************************************************');
fprintf (fid, '\n   PROGRAM INPUT FOR THE YAGI UDA ARRAY');
fprintf (fid, '\n   ******************************************************\n');

fprintf (fid, '\n   NUMBER OF MODES PER ELEMENT = %2d', M);
fprintf (fid, '\n   NUMBER OF ELEMENTS = %2d\n', N);

%     Print out the length of each element
a = 1;
while a <= (N-2) 
    fprintf (fid, '\n   LENGTH OF DIRECTOR # %2d = %12.5f WAVELENGTHS', a, L (a));
    a =  a + 1;
end
fprintf (fid, '\n   LENGTH OF REFLECTOR = %12.5f WAVELENGTHS', L (N-1));
fprintf (fid, '\n   LENGTH OF DRIVEN ELEMENT = %12.5f WAVELENGTHS\n', L (N));

%     Print out the element separation
fprintf (fid, '\n   SEPARATION BETWEEN DRIVEN ELEMENT & 1ST DIRECTOR = %12.5f WAVELENGTHS', S (1));
a = 2;
while a <= (N-2)
    fprintf (fid, '\n   SEPARATION BETWEEN DIRECTORS # %2d AND # %2d = %12.5f WAVELENGTHS', a-1, a, S (a));
    a =  a + 1;
end
fprintf (fid, '\n   SEPARATION BETWEEN REFLECTOR & DRIVEN ELEMENT = %12.5f WAVELENGTHS\n', S (N-1));

fprintf (fid, '\n   RADIUS FOR ALL ELEMENTS USED = %12.5e WAVELENGTHS\n\n\n', ALPHA);

%     Warn the user the program takes time to execute.
fprintf (1, '   WAIT ... This may take a few minutes !!!\n\n\n');

%  Open all necessary files for writing
fid1 = fopen ('Epl_yagi.dat', 'wt');
fid2 = fopen ('Hpl_yagi.dat', 'wt');
fid3 = fopen ('Cur_yagi.dat', 'wt');
fid4 = fopen ('Coe_yagi.dat', 'wt');


%     ************************************************************************************
%     ************************************************************************************
%     The real computation starts here.


%  Initialize some variables
a = 1;
while a <= (N - 2)
    YP (a) = a * S (a);
    a = a + 1;
end
YP (N-1) = - S (N-1);
YP (N) = 0;

RES = 0;
G2 = 0;
INDEX = 0;
DZ = L / (2 * M - 1);
ETA = 120 * pi;
MU = 4 * pi * 10 ^ (-7);
C = 3 * 10 ^ 8;
K = 2 * pi;
RTOD = 180 / pi;
DTOR = pi / 180;

A = zeros (M * N, M * N);
B = 1:(M*N);
B = B * 0;
Inm = zeros (N, M);

I = 1;
h=waitbar(0,'Program is running ...');
while I <= (M * N)
    waitbar(0.8*I/(M*N),h);
    IFACT = floor ((I - 1) / M);%     This determines the position of the observer with
    N1 = IFACT + 1;             %  N1 being the element at which the observer is.
    IMODE = I - IFACT * M;
    
    Z = (M - IMODE) * DZ (N1);  %  Based on the mode # and element find the distnce Z
    
    J = 1;
    while J <= (M * N)           
    
        JFACT = floor ((J - 1) / M); %     This determines the position of the source and its
        N2 = JFACT + 1;         %  corresponding mode.  N2 is the element on which the
        NMODE = J - JFACT * M;  %  source is and NMODE its mode #.
        
        if (N1 == N2)           %     If the effect of a mode is found on the element
            RHO = ALPHA;        %  that it is located, Y is the radius of the element.
        else                    %  Otherwise the distance Y is found using the formula
            RHO = YP (N1) - YP (N2);    %  (assuming X, X' = 0).
        end
        
        LL = 0;                 %     Define the limits of the integration
        UL = L (N2) / 2;

        RES = SINTEG (UL, LL, 10);  %     Perform numerical integration
        LEN = L (N2) / 2;
        G2 = KERNEL (LEN);
        F2M = NMODE * 2 - 1;

        A (I, J) = ETA / (j * 8 * pi ^ 2) * ...
                     ((F2M * pi / L (N2)) * (-1) ^ (NMODE + 1) * G2 + ...
                      (K ^ 2 - F2M ^ 2 * pi ^ 2 / L (N2) ^2) * RES);
        J = J + 1;

    end
    
    I = I + 1;
end

%  Fill the last row of the matrix corresponding to the feeder.
I = zeros (1, M * (N - 1));
J = ones (1, M);
A (M * N, :) = [I J];
B (M * N) = 1;

%  Invert the system to solve for the current coefficients in the
%  Fourier Series expansion.

ISIZE = N * M;
[A, IPERM, PIVOT] = LUDEC (A, ISIZE);
B = LUSOLV (A, ISIZE, IPERM, B);

%
% CONVERT THE SINGLE ARRAY OF THE CURRENT COEFFICIENTS TO A
% DOUBLE ARRAY OF THE FORM Imn.
%
NCUT = 0;

I = 1;
while I <= N
    
    J = 1;
    while J <= M
        Inm (I, J) = B (J + NCUT);
        J = J + 1;
    end
          
    NCUT = NCUT + M;
    I = I + 1;
end

%
% CALCULATE THE RADIATED FIELDS IN THE E-PLANE
%
% IN THIS PLANE THETA VARIES FROM 0 < THETA < 180, WHEREAS
% PHI = 90 IN HALF OF THE PATTERN AND PHI=270 IN THE OTHER 
% HALF.  THE PATTERN IS CALCULATED AT ONE DEGREE INCREMENTS. 
% 
NCUT = 0;

ML = 1;
% h=waitbar(0,'Program is running ...');
while ML <= 2
    
    if (ML == 1)
        PHI = 90 * DTOR;
        MAX = 181;
    else
        PHI = 270 * DTOR;
        MAX = 180;
    end
    
    ICOUNT = 1;
    while ICOUNT <= MAX
         waitbar(ICOUNT/MAX*ML*0.5*0.2+0.8,h);
        THETA = (ICOUNT - 1) * DTOR;
        
        if (THETA > pi)
            PHI = 270 * DTOR;
        end
        
        EZP = 0;
        
        I = 1;
        while I <= N
            
            IZP = 0;
            
            J = 1;
            while J <= M
                MODE = J;
                LEN = L (I);
                ANG = THETA;
                IZP = IZP + Inm (I, J) * ...
                            (ZMINUS (ANG, LEN, MODE) + ...
                             ZPLUS (ANG, LEN, MODE));
                J = J + 1;
            end
            
            AEXP = K * YP (I) * sin (THETA) * sin(PHI);
            EZP = EZP + L (I) * exp (j * AEXP) * IZP;
            
            I = I + 1;
        end
        
        ETHETA (NCUT + ICOUNT) = j * C * MU / 8 * sin (THETA) * EZP;

        ICOUNT = ICOUNT + 1;
    end
    
    NCUT = NCUT + MAX;
    ML = ML + 1;
end
close(h);

%  FIND THE MAXIMUM VALUE IN THE E-PLANE PATTERN
EMAX = 10 ^ (-12);

abs_ETHETA = abs (ETHETA);
ARG = max (abs_ETHETA);
if ARG > EMAX
    EMAX = ARG;
end

% NORMALIZE THE PATTERN TO THE MAXIMUM VALUE, CONVERT IN dB, 
% AND WRITE OUT THE RESULTS

fprintf (fid1, '# E-PLANE PATTERN OF THE YAGI UDA ANTENNA\n');
fprintf (fid1, '# =======================================\n#\n');
fprintf (fid1, '#      THETA      E-THETA (THETA, PHI=90 OR 270)\n#\n');

I = 1;
while I <= 361
    THETA = I - 1;
    ARG = abs (ETHETA (I));
    
    if ((ARG/EMAX) > (10 ^ (-6)))
        ETH (I) = 20 * log10 (ARG / EMAX);
    else
        ETH (I) = -120;
    end
    
    fprintf (fid1, ' %12.4f%12.4f\n', THETA, ETH (I));
    I = I + 1;
end

% RECORD THE E-PLANE VALUES FOR FUTURE PLOTTING.
E_PLANE = ETH;

%
% FIND THE FRONT-TO-BACK RATIO IN THE E-PLANE PATTERN 
%
EFTOB = - ETH (271);

%
% FIND THE 3-dB BEAMWIDTH IN THE E-PLANE PATTERN
%
I = 91;
while I <= 270
    ETH (I) = ETH (I) + 3;
    I = I + 1;
end

I = 91;
while I <= 270
    THETA = I-1;
    
    if (ETH (I) == 0)
        E3D_BW = 2 * ((I - 1) - 90);
        I = 300;    % Terminate while loop early
    elseif ((ETH (I - 1) > 0) & (ETH (I) < 0))
        E3D_BW = 2 * (- ETH (I) / (ETH (I) - ETH (I - 1)) + (I - 1) - 90);
        I = 300;
    end
    
    I = I + 1;
end

        
%
% CALCULATE THE RADIATED FIELDS IN THE H-PLANE 
%
% IN THIS PLANE THETA = 90 AND PHI VARIES FROM 0 < PHI < 360.
%
THETA = 90 * DTOR;
MAX = 361;

ICOUNT = 1;
while ICOUNT <= MAX
    
    PHI = (ICOUNT - 1) * DTOR;
    EZP = 0;
    
    I = 1;
    while I <= N
        IZP = 0;
        
        J = 1;
        while J <= M
            MODE = J;
            LEN = L (I);
            ANG = PHI;
            IZP = IZP + Inm (I, J) * ...
                        (ZMINUS (ANG, LEN, MODE) + ...
                         ZPLUS (ANG, LEN, MODE));
            J = J + 1;
        end
        
        AEXP = K * YP (I) * sin (THETA) * sin(PHI);
        EZP = EZP + L (I) * exp (j * AEXP) * IZP;
        I = I + 1;
    end
    
    ETHETA (ICOUNT) = j * C * MU / 8 * sin (THETA) * EZP;
    ICOUNT = ICOUNT + 1;
end

%
%  FIND THE MAXIMUM VALUE IN THE H-PLANE PATTERN
%
EMAX = 10 ^ (-12);
abs_ETHETA = abs (ETHETA);
ARG = max (abs_ETHETA);
if (ARG > EMAX)
    EMAX = ARG;
end

%
% NORMALIZE THE PATTERN TO THE MAXIMUM VALUE, CONVERT IN dB, AND 
% WRITE OUT THE RESULTS
%
fprintf (fid2, '# H-PLANE PATTERN OF THE YAGI UDA ANTENNA\n');
fprintf (fid2, '# =======================================\n#\n');
fprintf (fid2, '#       PHI       E-THETA (PHI, THETA=90)\n#\n');

I = 1;
while I <= 361
    
    PHI = I - 1;
    ARG = abs (ETHETA (I));
    
    if (ARG / EMAX) > (10 ^ (-6))
        ETH (I) = 20 * log10 (ARG / EMAX);
    else
        ETH (I) = - 120;
    end
    
    fprintf (fid2, ' %12.4f%12.4f\n', PHI, ETH (I));
    I = I + 1;
end

% RECORD THE H-PLANE VALUES FOR FUTURE PLOTTING
H_PLANE = ETH;

%
% FIND THE FRONT-TO-BACK RATIO IN THE H-PLANE PATTERN 
%
HFTOB = - ETH (271);

%
% FIND THE 3-dB BEAMWIDTH IN THE H-PLANE PATTERN
%
I = 1;
while I <= 181
    ETH (I) = ETH (I) + 3;
    I = I + 1;
end

I = 91;
while I <= 270
    PHI = I - 1;
    
    if (ETH (I) == 0)
        H3D_BW = 2 * ((I - 1) - 90);
        I = 300;
    elseif (ETH (I-1) > 0) & (ETH (I) < 0)
        H3D_BW = 2 * (- ETH (I) / (ETH (I) - ETH (I - 1)) + (I - 1) - 90);
        I = 300;
    end
    
    I = I + 1;
end

%
% CALCULATE THE ANTENNA DIRECTIVITY
%
THETA = 90 * DTOR;
PHI = 90 * DTOR;
AZ = 0;

I = 1;
while I <= N
    IZP = 0;
    
    J = 1;
    while J <= M
        MODE = J;
        LEN = L (I);
        ANG = THETA;
        IZP = IZP + Inm (I, J) * ...
                    (ZMINUS (ANG, LEN, MODE) + ...
                     ZPLUS (ANG, LEN, MODE));
        J = J + 1;
    end
    
    AEXP = K * YP (I) * sin (THETA) * sin (PHI);
    AZ = AZ + L (I) * exp (j * AEXP) * IZP;
    
    I = I + 1;
end

UMAX = 3.75 * pi * abs (AZ) ^ 2 * sin (THETA) ^ 2;

PRAD = SCINT2 (0, pi, 0, 2 * pi, N, M, Inm, YP);

D0 = 4 * pi * UMAX / abs (PRAD);

fprintf (fid, '\n   ******************************************************\n');
fprintf (fid, '   PROGRAM OUTPUT FOR THE YAGI UDA ARRAY\n');
fprintf (fid, '   ******************************************************\n');
fprintf (fid, '\n   3-dB BEAMWIDTH IN THE E-PLANE PATTERN = %12.2f  DEGREES\n', E3D_BW);
fprintf (fid, '\n   3-dB BEAMWIDTH IN THE H-PLANE PATTERN = %12.2f  DEGREES\n', H3D_BW);
fprintf (fid, '\n   FRONT-TO-BACK RATIO IN THE E-PLANE = %12.4f  dB\n', EFTOB);
fprintf (fid, '\n   FRONT-TO-BACK RATIO IN THE H-PLANE = %12.4f  dB\n', HFTOB);
fprintf (fid, '\n   DIRECTIVITY = %16.3f  dB\n', 10 * log10 (D0));

%     BASED ON THE FOURIER COEFFICIENTS OF THE CURRENT, CALCULATE THE 
%     CURRENT DISTRIBUTION ON THE ELEMENTS.  NOTE THAT EACH ELEMENT 
%     IS SUBDIVIDED INTO 100 SECTIONS FOR THIS CALCULATION.
%
%
fprintf (fid3, '   CURRENT DISTRIBUTION ON THE ELEMENTS\n');
fprintf (fid3, '   ====================================\n\n');

IL = 1;
while IL <= N
    
    fprintf (fid3, '\n\n   ELEMENT # %3d\n', IL);
    fprintf (fid3, '   =========\n\n');
    
    DZ (IL) = L (IL) / 100;
    
    I = 1;
    while I <= 51
        
        Z = (I - 1) * DZ (IL);
        IZP = 0;
        
        J = 1;
        while J <= M
            
            F2M = 2 * J - 1;
            IZP = IZP + Inm (IL, J) * cos (F2M * pi * Z / L (IL));
            
            J = J + 1;
        end
        
        CUR (I) = abs (IZP);
        angle = atan2 (imag (IZP), real (IZP)); 
        PHA (I) = angle * RTOD;
        
        I = I + 1;
    end

    fprintf (fid3, '       DISTANCE        CURRENT     CURRENT\n');
    fprintf (fid3, '                      MAGNITUDE     PHASE \n\n');

    I = 1;
    while I <= 51
        Z = (I - 1) * DZ (IL);
        fprintf (fid3, '   %12.5f   %12.5f   %12.5f\n', Z, CUR (I), PHA (I));
        I = I + 1;
    end

    CENTER_CURRENT (IL) = CUR (1); % Record for future plotting
    IL = IL + 1;
end

I = 1;
while I <= N
    
    fprintf (fid4, '\n\n   ELEMENT # %3d\n', I);
    fprintf (fid4, '   =========\n\n');
    fprintf (fid4, '  MODE #    MAGNITUDE      PHASE\n\n');
    
    J = 1;
    while J <= M
        
        CURRENT = abs (Inm (I, J));
        angle_radian = atan2 (imag (Inm (I, J)), real (Inm (I, J)));
        ANGLE = angle_radian * RTOD;
        
        fprintf (fid4, '   %2d   %12.5f   %12.5f\n', J, CURRENT, ANGLE);
        
        J = J + 1;
    end
    
    I = I + 1;
end

I = 1;
while I > 0
    
    fprintf (fid, '\n\n   *** NOTE:\n');
    fprintf (fid, '       E-PLANE PATTERN IS STORED IN Epl-yagi.dat\n');
    fprintf (fid, '       H-PLANE PATTERN IS STORED IN Hpl-yagi.dat\n');
    fprintf (fid, '       CURRENT ON EACH ELEMENT IS STORED IN Cur-yagi.dat\n');
    fprintf (fid, '       MODE COEFFs. FOR EACH ELEMENT ARE STORED IN Coe-yagi.dat\n\n');

    if fid == 1
        I = 0;  % exit the loop
    else
        fid = 1;    % print same message on the screen again
    end
    
end

fclose ('all');

E_PLANE = E_PLANE (1:360);
H_PLANE = H_PLANE (1:360);
angle = 1:1:360;

figure;
plot (angle, E_PLANE, '-b', 'LineWidth', 2);
hold on;
plot (angle, H_PLANE, '--r', 'LineWidth', 2);
legend ('E-Plane', 'H-Plane');
xlim ([1 360]);
ylim ([-60 0]);
title ('Yagi-Uda Analysis');
xlabel ('Theta(E)/Phi(H) degrees');
ylabel ('Field Pattern (dB)');
hold off;

figure;
INDEX = 1:1:N;
if N >= 3
    CENTER_CURRENT = [CENTER_CURRENT(N-1:N) CENTER_CURRENT(1:N-2)];
    INDEX = [INDEX(N-1:N) INDEX(1:N-2)]; 
end
plot (CENTER_CURRENT, 'LineWidth', 2);
set (gca, 'XTick', 1:1:N);
set (gca, 'XTickLabel',INDEX);
title ('Current Distribution');
ylim ([0 1]);
xlabel ('Element Number');
ylabel ('Element Current Amplitude');

figure;
h1=elevation(angle*pi/180,E_PLANE,-40,0,5,'b'); hold on;
h2=elevation(angle*pi/180,H_PLANE,-40,0,5,'r--');
set([h1 h2],'linewidth',2);
legend([h1 h2],{'E-Plane','H-Plane'});
    
%  End of the yagi_uda function
%***************************************************************** 


%***************************************************************** 
%     FUNCTION SINTEG  (SINGLE PRECISION)
%
%     PURPOSE
%     TO PERFORM COMPLEX SINGLE INTEGRATION
%     DOES SIXTEEN POINT GAUSSIAN QUADRATURE INTEGRATION
%     WITH INCREASING ACCURACY SET BY INTEGER NO
%
%     USAGE
%     ANS = SINTEG(UL,LL,NO)
%
%     DESCRIPTION OF PARAMETERS
%     UL  -  UPPER LIMIT OF THE INTEGRATION (REAL)
%     LL  -  LOWER LIMIT OF THE INTEGRATION (REAL)
%     NO  -  NUMBER OF DIVISIONS BETWEEN LL AND UL  (INTEGER)
%
%     ANS -  INTEGRATION RESULT
%     L   -  THE LENGTHS OF THE ARRAY ELEMENTS (GLOBAL VARIABLE)

function ANS = SINTEG (UL, LL, NO)

GAUSS = [-0.0950125098376370
         -0.2816035507792590
         -0.4580167776572270
         -0.6178762444026440
         -0.7554044083550030
         -0.8656312023878320
         -0.9445750230732330
         -0.9894009349916500
         0.0950125098376370
         0.2816035507792590
         0.4580167776572270
         0.6178762444026440
         0.7554044083550030
         0.8656312023878320
         0.9445750230732330
         0.9894009349916500];
     
LEGEND = [0.1894506104550680
          0.1826034150449240
          0.1691565193950020
          0.1495959888165770
          0.1246289712555340
          0.0951585116824930
          0.0622535239386480
          0.0271524594117540
          0.1894506104550680
          0.1826034150449240
          0.1691565193950020
          0.1495959888165770
          0.1246289712555340
          0.0951585116824930
          0.0622535239386480
          0.0271524594117540];

DEL = (UL - LL) / (2 * NO);
SUM = 0;

J = 1;
while J <= NO
    S = LL + (2 * J - 1) * DEL;
    I = 1;
    while I <= 16
        X = S + GAUSS (I) * DEL;
        SUM = SUM + LEGEND (I) * FF (X);
        I = I + 1;
    end
    J = J + 1;
end

ANS = SUM * DEL;

%  End of the SINTEG function
%***************************************************************** 
     
      

%***************************************************************** 
%     FUNCTION FF
%
function value = FF (X)
global Z RHO N2 NMODE L

K = 2 * pi;

RMINUS = sqrt (RHO ^ 2 + (Z - X) ^ 2);
RPLUS = sqrt ( RHO ^ 2 + (Z + X) ^ 2);

F2M = NMODE * 2 - 1;

value = (exp (-j * K * RMINUS) / (4 * pi * RMINUS) + ...
         exp (-j * K * RPLUS) / (4 * pi * RPLUS)) * ... 
        cos (F2M * pi * X / L (N2));

%  End of the FF function
%***************************************************************** 


%****************************************************************
%     FUNCTION KERNEL
%
function G2 = KERNEL (ZP)
global Z RHO
     
K = 2 * pi;

RMINUS = sqrt (RHO ^ 2 + (Z - ZP) ^ 2);
RPLUS = sqrt (RHO ^ 2 + (Z + ZP) ^ 2);

G2 = exp (-j * K * RMINUS) / (4 * pi * RMINUS) + ...
     exp (-j * K * RPLUS) / (4 * pi * RPLUS);

%  End of the KERNEL function
%****************************************************************


%****************************************************************
%     FUNCTION ZMINUS
%
function value = ZMINUS (TH, LG, NM)

K = 2 * pi;

F2M = 2 * NM - 1;

ARG1 = (F2M * pi / LG - K * cos (TH)) * (LG / 2);
if (ARG1 == 0)
    value = 1;
else
    value = sin (ARG1) / ARG1;
end

%  End of the ZMINUS function
%****************************************************************


%****************************************************************
%     FUNCTION ZPLUS
%
function value = ZPLUS (TH, LG, NM)

K = 2 * pi;

F2M = 2 * NM - 1;

ARG2 = (F2M * pi / LG + K * cos (TH)) * (LG / 2);
if (ARG2 == 0)
    value = 1;
else
    value = sin (ARG2) / ARG2;
end

%  End of the ZPLUS function
%****************************************************************


%****************************************************************
%     FUNCTION FXZ
%
function value = FXZ (THETA, PHI, N, M, Inm, YP)
global L

K = 2 * pi;
AZ = 0;

I = 1;
while I <= N
    
    IZP = 0;
    
    J = 1;
    while J <= M
        
        MODE = J;
        LEN = L (I);
        ANG = THETA;
        
        IZP = IZP + Inm (I, J) * ...
                    (ZMINUS (ANG, LEN, MODE) + ...
                     ZPLUS (ANG, LEN, MODE));
        J = J + 1;         
    end
    
    AEXP = K * YP (I) * sin (THETA) * sin (PHI);
    AZ = AZ + L(I) * exp (j * AEXP) * IZP;
    
    I = I + 1;
end

U = 3.75 * pi * abs (AZ) ^ 2 * sin (THETA) ^ 2;

value = sin (THETA) * U;

%  End of the FXZ function
%****************************************************************



%****************************************************************
%     FUNCTION LUDEC
%
function [A, IPERM, SCAL] = LUDEC (Z, N)

%    REPLACES MATRIX BY ITS LU DECOMPOSITION
%    GET SCALING INFO.

I = 1;
while I <= N
    
    ZMAX = 0;
    
    J = 1;
    while J <= N
        
        CAZ = abs (Z (I,J));
        if (CAZ > ZMAX) 
            ZMAX = CAZ;
        end
        J = J + 1;
    end
    
    SCAL (I) = 1 / ZMAX;
    I = I + 1;
end

%    CROUT's algorithm.
J = 1;
while J <= N

    I = 1;
    while I <= (J - 1)
        
        K = 1;
        while K <= (I - 1)
            Z (I, J) = Z (I, J) - Z (I, K) * Z (K, J);
            K = K + 1;
        end
        
        I = I + 1;
    end

%    SEARCH FOR LARGEST PIVOT ELEMENT.
    ZMAX = 0;
    
    I = J;
    while I <= N
        
       K = 1;
       while K <= (J - 1)
           Z (I, J) = Z (I, J) - Z (I, K) * Z (K, J);
           K = K + 1;
       end
        
       TEST = SCAL (I) * abs (Z (I, J));
       
       if (round ((TEST - ZMAX) * 10 ^ 8) > 0)  % Work around way
           IMAX = I;
           ZMAX = TEST;
       end

       I = I + 1;
    end

%    INTERCHANGE THE ROWS.
    if (J ~= IMAX)
         
        K = 1;
        while K <= N
            TEMP = Z (IMAX, K);
            Z (IMAX, K) = Z (J, K);
            Z (J, K) = TEMP;
            K = K + 1;
        end         
        SCAL (IMAX) = SCAL (J);
    end
   
%   DIVIDE BY PIVOT ELEMENT.
    IPERM (J) = IMAX;
    if (J ~= N)
         
        I = J + 1;
        while I <= N
            Z (I, J) = Z (I, J) / Z (J, J);
            I = I + 1;
        end
    end
     
    J = J + 1;
end

A = Z;

%  End of the LUDEC function
%****************************************************************



%****************************************************************
%     FUNCTION LUSOLV
%
function B = LUSOLV (Z, N, IPERM, V)

%    SOLVES LINEAR SYSTEM GIVEN THE LU DECOMPOSITION FROM LUDEC
%    FORCING VECTOR IS REPLACED WITH SOLUTION VECTOR UPON EXIT 

%    FORWARD SUBSTITUTION.
I = 1;
while I <= N
    
    TEMP = V (IPERM (I));
    V (IPERM (I)) = V (I);
    
    J = 1;
    while J <= (I - 1)
        TEMP = TEMP - Z (I, J) * V (J);
        J = J + 1;
    end
    
    V (I) = TEMP;
    I = I + 1;
end

%    BACKWARD SUBSTITUTION.
I = 1;
while I <= N
    II = N - I + 1;
    TEMP = V (II);
    
    J = II + 1;
    while J <= N
        TEMP = TEMP - Z (II, J) * V (J);
        J = J + 1;
    end
    
    V (II) = TEMP / Z (II, II);
    I = I + 1;
end

B = V;
%  End of the LUSOLV function
%****************************************************************


%******************************************************************
% FUNCTION SCINT2
% SCINT2 IS A SINGLE PRECISION, COMPLEX, INTEGRATION ROUTINE
% IN 2 DIMENSIONS.  THIS ROUTINE USES 16 POINT GAUSSIAN QUADRATURES,
% WITH LEGANDRE COEFFICIENTS.  ENTER WITH:
%
%       (R)X1          LOWER LIMIT OF X INTEGRATION
%       (R)X2          UPPER LIMIT OF X INTEGRATION
%       (R)Z1          LOWER LIMIT OF Z INTEGRATION
%       (R)Z2          UPPER LIMIT OF Z INTEGRATION
%       (C)RES         RESULTS OF INTEGRATION
%       (I)N           NUMBER OF ELEMENTS IN THE ARRAY
%       (I)M           NUMBER OF MODES PER ELEMENT
%       (R)Inm         THE ARRAY OF CURRENT COEFFICIENTS
%       (R)YP          THE ELEMENT DISTANCES ALONG THE Y AXIS
%
%***********************************************************************
%
%
function RES = SCINT2 (X1, X2, Z1, Z2, N, M, Inm, YP)

R = [0.0950125098
     0.2816035508
     0.4580167777
     0.6178762444
     0.7554044084
     0.8656312024
     0.9445750231
     0.9894009350
     -0.9894009350
     -0.9445750231
     -0.8656312024
     -0.7554044084
     -0.6178762444
     -0.4580167777
     -0.2816035508
     -0.0950125098];
 
W = [0.1894506105
     0.1826034150
     0.1691565194
     0.1495959888
     0.1246289713
     0.0951585117
     0.0622535239
     0.0271524594
     0.0271524594
     0.0622535239
     0.0951585117
     0.1246289713
     0.1495959888
     0.1691565194
     0.1826034150
     0.1894506105];
 
SX = 0.5 * (X2 + X1);
DX = 0.5 * (X2 - X1);
SZ = 0.5 * (Z2 + Z1);
DZ = 0.5 * (Z2 - Z1);

TT = 0;

J = 1;
while J <= 16
    
    Z  = SZ + DZ * R (J);
    SS = 0;
    
    I = 1;
    while I <= 16
        
        X  = SX + DX * R (I);
        SS = W (I) * FXZ (X, Z, N, M, Inm, YP) + SS;
        
        I = I + 1;
    end
    
    S = DX * SS;
    TT = S * W (J) + TT;
    
    J = J + 1;
end

RES = DZ * TT;
%  End of the SINT2 function
%****************************************************************
