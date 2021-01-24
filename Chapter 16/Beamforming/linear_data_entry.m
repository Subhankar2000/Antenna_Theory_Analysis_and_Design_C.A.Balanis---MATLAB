%************************************************************************
%	[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry
%*************************************************************************
%	DATA_ENTRY is a MATLAB function that collects the user's input data in
%	question form.
%
%	Output Parameters Description
%	-----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%   - sig        signal parameters (amplitudes and directions)
%                for SOI and SNOI
%   - noise      noise parameter (mean and variance)
%   - type       sinusoidal or BPSK signal
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - AF_thresh  nulls' depth
%   - Mu         convergence factor for LMS algorithm
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function [N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry

%%%%%%%%%%%%%%%%%%% Default Values %%%%%%%%%%%%%%%%%%
def_N              = 8;
def_d              = 0.5;
def_SOI            = 1;
def_q              = 1;
def_SNOI           = 1;
def_noise_mean     = 0;
def_noise_var      = 0.1;
def_type_n         = 2;
def_AF_thresh      = -40;
def_Mu             = 0.001;
def_nn             = 500;
def_E_pattern_file = 'linear_isotropic.e';
%---------------------------------------------------%

%%%%%%%%%%%%% Strings initialization %%%%%%%%%%%%%%%%
N_string           = sprintf('Enter number of elements in linear smart antenna [%d]: ',def_N);
d_string           = sprintf('Enter the spacing d (in lambda) between adjacent elements [%2.1f]: ',def_d);
SOI_string         = sprintf('Enter the Pilot signal (SOI) amplitude [%d]: ',def_SOI);
SOId_string        = sprintf('Enter the Pilot signal (SOI) direction (degrees between 0 and 90): ');
q_string           = sprintf('Enter number of interfering signals (SNOI) [%d]: ',def_q);
noise_mean_string  = sprintf('Enter the mean of the noise [%d]: ',def_noise_mean);
noise_var_string   = sprintf('Enter the variance of noise [%2.1f]: ',def_noise_var);
type               = sprintf('Type of signal:\n\t[1] sinusoid\n\t[2] BPSK\nEnter number [%d]: ',def_type_n);
nn_string          = sprintf('Enter the number of data samples [%d]: ',def_nn); 
AF_thresh_string   = sprintf('Enter AF threshold (dB) [%d]: ',def_AF_thresh);
Mu_string          = sprintf('Enter a value for Mu of LMS algorithm (0 < Mu < 1) [%5.4f]: ',def_Mu);
E_pattern_string   = sprintf('Enter element pattern filename (*.e) [%s]: ',def_E_pattern_file);
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%% Error Messages %%%%%%%%%%%%%%%%%%%%
err_1             = sprintf('\nSignal type not supported...');
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%
% ------------ Number of elements? ------------- %
N = input(N_string);
if isempty(N)
   N = def_N;
end;
% ---------------------------------------------- %

% ------------ Inter-element spacing? ---------- %
d = input(d_string);
if isempty(d)
   d = def_d;
end;
% ---------------------------------------------- %
    
% -------------- SOI amplitude? ---------------- %
SOI = input(SOI_string);
if isempty(SOI)
   SOI = def_SOI;
end;
% ---------------------------------------------- %
    
% -------------- SOI direction? ---------------- %
SOId = [];
while isempty(SOId)
    SOId = input(SOId_string);
end;
% ---------------------------------------------- %
    
% ------ SNOIs amplitudes and directions? ------ %
q = input(q_string);
if isempty(q)
   q = def_q;
end;
    
sig = [SOI SOId];
    
for k = 1 : q,
    SNOI_k_string      = sprintf('Enter the amplitude of No. %d Interference signal (SNOI_%d) [%d]: ',k,k,def_SNOI);
    SNOId_k_string     = sprintf('Enter the direction of No. %d Interference signal (SNOI_%d) (degrees between 0 and 90): ',k,k);
    SNOI_k = input(SNOI_k_string);
    if isempty(SNOI_k)
       SNOI_k = def_SNOI;
    end;
    SNOId_k = [];
    while isempty(SNOId_k)
        SNOId_k = input(SNOId_k_string);
    end;
    sig = [sig; SNOI_k SNOId_k];
end;
% ---------------------------------------------- %
    
% ----------------- Noise data? ---------------- %
noise_string = input('Insert noise? ([y]/n):','s');
if (isempty(noise_string) | noise_string == 'y')
   noise_mean = input(noise_mean_string);
   if isempty(noise_mean)
      noise_mean = def_noise_mean;
   end;
   noise_var = input(noise_var_string);
   if isempty(noise_var)
      noise_var = def_noise_var;
   end;
   noise = [noise_mean noise_var];
else
   fprintf('------------No noise is inserted.-----------\n');
   noise = [];
end;
% ---------------------------------------------- %
    
% ----------------- Signal type? --------------- %
type_n = [];
if isempty(type_n)
   type_n = def_type_n;
end;
switch type_n
case 1
    type = 'sinusoid';
    def_NN = 100;
case 2
    type = 'bpsk';
    def_NN = 1;
otherwise
    error(err_1);
end;
% ---------------------------------------------- %
    
% ------------ Number of data samples? --------- %
nn = input(nn_string);
if isempty(nn)
    nn = def_nn;
end;
% ---------------------------------------------- %

% ------- Number of samples per symbol? -------- %
NN = [];
if isempty(NN)
   NN = def_NN;
end;
% ---------------------------------------------- %

% ---------------- Mu for LMS? ----------------- %
Mu = input(Mu_string);
if isempty(Mu)
   Mu = def_Mu;
end;
% ---------------------------------------------- %

% ---------------- Nulls depth? ---------------- %
if q==0
   AF_thresh = def_AF_thresh;
else
   AF_thresh = input(AF_thresh_string);
   if isempty(AF_thresh)
      AF_thresh = def_AF_thresh;
   end;
end;
% ---------------------------------------------- %

% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
   E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;
