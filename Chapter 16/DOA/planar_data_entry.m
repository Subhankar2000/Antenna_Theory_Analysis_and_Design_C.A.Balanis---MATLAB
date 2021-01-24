%************************************************************************
%	[N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry
%***********************************************************************
%	DATA_ENTRY is a MATLAB function that collects the
%   user's input data in question form
%
%	Output Parameters Description
%	-----------------------------
%	- N          number of elements in x axis and y axis
%	- d          inter-element spacing (in wavelength) in x and y axis
%                default value is 0.5 (Nyquist rate)
%   - sig        signal parameters (amplitudes and directions)
%                for SOI and SNOI
%   - noise      noise parameter (mean and variance)
%   - type       sinusoidal or BPSK signal
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function [N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry

%%%%%%%%%%%%%%%%%%% Default Values %%%%%%%%%%%%%%%%%%
def_N                   = [8 8];
def_d                   = [0.5 0.5];
def_SOI                 = 1;
def_q                   = 1;
def_SNOI                = 1;
def_noise_mean          = 0;
def_noise_var           = 0.1;
def_type_n              = 2;
def_nn                  = 500;
def_E_pattern_file      = 'planar_isotropic.e';
%---------------------------------------------------%
   
%%%%%%%%%%%%% Strings initialization %%%%%%%%%%%%%%%%
N_string_x              = sprintf('Enter number of elements in the x-axis [%d]: ',def_N(1));
N_string_y              = sprintf('Enter number of elements in the y-axis [%d]: ',def_N(2));
d_string_x              = sprintf('Enter the spacing in the x-axis, dx (lambda) [%2.1f]: ',def_d(1));
d_string_y              = sprintf('Enter the spacing in the y-axis, dy (lambda) [%2.1f]: ',def_d(2));
SOI_string              = sprintf('Enter the Pilot signal (SOI) amplitude [%d]: ',def_SOI);
SOId_string_theta       = sprintf('Enter the Pilot signal (SOI) in theta direction (degrees between 0 and 90): ');
SOId_string_phi         = sprintf('Enter the Pilot singal (SOI) in phi direction (degrees between 0 and 360): ');
q_string                = sprintf('Enter number of interfering signals (SNOI) [%d]: ',def_q);
noise_mean_string       = sprintf('Enter the mean of the noise [%d]: ',def_noise_mean);
noise_var_string        = sprintf('Enter the variance of the noise [%2.1f]: ',def_noise_var);
type                    = sprintf('Type of signal:\n\t[1] sinusoid\n\t[2] BPSK\nEnter number [%d]: ',def_type_n);
nn_string               = sprintf('Enter the number of data samples [%d]: ',def_nn); 
E_pattern_string        = sprintf('Enter element pattern filename (*.e) [%s]: ',def_E_pattern_file);
%---------------------------------------------------%

%%%%%%%%%%%%%%%%% Error Messages %%%%%%%%%%%%%%%%%%%%
err_1                   = sprintf('\nSignal type is not supported...');
err_2                   = sprintf('\nAlgorithm type is not supported...');
err_3                   = sprintf('\nMethod is not supported...');
err_4                   = sprintf('\nThis particular weights initialization Method is not supported...');
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%
% ------------ Number of elements? ------------- %
N_x = input(N_string_x);
if isempty(N_x)
    N_x = def_N(1);
end;
N_y = input(N_string_y);
if isempty(N_y)
    N_y = def_N(2);
end;
N = [N_x N_y];
% ---------------------------------------------- %

% ------------ Inter-element spacing? ---------- %
d_x = input(d_string_x);
if isempty(d_x)
    d_x = def_d(1);
end;
d_y = input(d_string_y);
if isempty(d_y)
    d_y = def_d(2);
end;
d = [d_x d_y];
% ---------------------------------------------- %
    
% -------------- SOI amplitude? ---------------- %
SOI = input(SOI_string);
if isempty(SOI)
    SOI = def_SOI;
end;
% ---------------------------------------------- %
    
% -------------- SOI direction? ---------------- %
SOId_theta = []; SOId_phi = [];
while isempty(SOId_theta)
    SOId_theta = input(SOId_string_theta);
end;
while isempty(SOId_phi)
    SOId_phi = input(SOId_string_phi);
end;
% ---------------------------------------------- %
    
% ------ SNOIs amplitudes and directions? ------ %
q = input(q_string);
if isempty(q)
    q = def_q;
end;
    
sig = [SOI SOId_theta SOId_phi];
    
for k = 1 : q,
    SNOI_k_string           = sprintf('Enter the amplitude of No. %d Interference singal (SNOI_%d)  [%d]: ',k,k, def_SNOI);
    SNOId_k_string_theta    = sprintf('Enter the theta direction of No. %d Interference singal (SNOI_%d) (degrees between 0 and 90): ',k,k);
    SNOId_k_string_phi      = sprintf('Enter the phi direction of No. %d Interference singal (SNOI_%d) (degrees between 0 and 360): ',k,k);
    SNOI_k = input(SNOI_k_string);
    if isempty(SNOI_k)
        SNOI_k = def_SNOI;
    end;
    SNOId_k_theta = []; SNOId_k_phi = [];
    while isempty(SNOId_k_theta)
        SNOId_k_theta = input(SNOId_k_string_theta);
    end;
    while isempty(SNOId_k_phi)
        SNOId_k_phi = input(SNOId_k_string_phi);
    end;
    sig = [sig; SNOI_k SNOId_k_theta SNOId_k_phi];
end;
% ---------------------------------------------- %
    
% ----------------- Noise data? ---------------- %
noise_string = input('Insert noise? ([y]/n): ','s');
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

% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
    E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;







