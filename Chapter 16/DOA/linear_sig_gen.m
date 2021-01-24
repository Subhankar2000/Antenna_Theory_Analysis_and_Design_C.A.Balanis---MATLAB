%**********************************************************************
%	DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)
%**********************************************************************
%	SIG_GEN is a MATLAB function that generates the SOI and SNOIs
%   signals (with and without Gaussian noise)
%
%	Input Parameters Description
%	----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - type       option: 'sinusoid' or 'bpsk'
%   - sig        signals amplitudes and directions of SOI and SNOI
%   - noise      amplitude and variance values
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%
%	Output Parameters Description
%	-----------------------------
%   -DOAs        Estimations of DOAs
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

function DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)

%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
k0 = 2*pi;                     % k in free space
n_sig = size(sig,1);           % Total number of signals
m = [];

x = linspace(0,N-1,N);

y = 0;
sig(:,1) = sig(:,1) .* E_pattern(length(E_pattern) - ((length(E_pattern) - 1) / 2 - sig(:,2))).';
sigr = [sig(:,1) pi/180*sig(:,2) zeros(n_sig,1)];       % degrees to radians conversion

[X,Y] = meshgrid(x,y);

err_1 = sprintf('\nSignal type not supported...');
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SOI %%%%%%%%%%%%%%%%%%%%
PSI_SOI = -k0*d*sin(sigr(1,2))*cos(sigr(1,3))*X - k0*d*sin(sigr(1,2))*sin(sigr(1,3))*Y;

switch type
   case 'sinusoid'
      t = linspace(0,nn-1/NN,NN*nn);
      [T1,PSI_SOI] = meshgrid(t,PSI_SOI);
      s = sigr(1,1)*exp(i*(2*pi*T1+PSI_SOI)); 
      S_SOI  = real(s(1,:));
   case 'bpsk'
      STATE1 = sum(100*clock);
      rand('state',STATE1);
      s_rand = round(rand(nn,1));
      s_rand_index = find(s_rand == 0);
      s_rand(s_rand_index) = -1;
      s_NN = ones(1,NN);
      s_rand = (s_rand*s_NN)';
      s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
      [s_rand,PSI_SOI] = meshgrid(s_rand,PSI_SOI);
      s = sigr(1,1)*(s_rand.*exp(i*PSI_SOI));
      S_SOI  = s(1,:);
   otherwise
      error(err_1);
end;

if ~isempty(noise)
    s_noise = s;

    for k = 1 : N
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s_noise(k,:) = s_noise(k,:) + noise_data;
    end;
    DOA_theta =est_doa(s_noise,N,1,d,1,1);
else
    DOA_theta =est_doa(s,N,1,d,1,1);
end;    
DOAs = DOA_theta;
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SNOI %%%%%%%%%%%%%%%%%%%
for k = 2 : n_sig,
   PSI_SNOI = -k0*d*sin(sigr(k,2))*cos(sigr(k,3))*X - k0*d*sin(sigr(k,2))*sin(sigr(k,3))*Y;
   switch type
      case 'sinusoid'
         t = linspace(0,nn-1/NN,NN*nn);
         [T2,PSI_SNOI] = meshgrid(t,PSI_SNOI);
         m = [1; 100*rand(n_sig-1,1)];                        % frequency multipliers
         SNOI = sigr(k,1)*exp(i*(2*pi*m(k)*T2+PSI_SNOI));
         s = s + SNOI;
      case 'bpsk'
         STATE2 = sum(k*100*clock);
         rand('state',STATE2);
         s_rand = round(rand(nn,1));
         s_rand_index = find(s_rand == 0);
         s_rand(s_rand_index) = -1;
         s_NN = ones(1,NN);
         s_rand = (s_rand*s_NN)';
         s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
         [s_rand,PSI_SNOI] = meshgrid(s_rand,PSI_SNOI);
         SNOI = sigr(k,1)*(s_rand.*exp(i*PSI_SNOI));
         if ~isempty(noise)
             s_noise = SNOI;

             for k = 1 : N
                 STATE3 = sum(rand(1)*100*clock);
                 randn('state',STATE3);
                 noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 STATE4 = sum(rand(1)*100*clock);
                 randn('state',STATE4);
                 noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 noise_data = complex(noise_data_real,noise_data_imag);
                 s_noise(k,:) = s_noise(k,:) + noise_data;
             end;
             DOA_theta = est_doa(s_noise,N,1,d,1,1);
         else
             DOA_theta = est_doa(SNOI,N,1,d,1,1);
         end;    
         DOAs = [DOAs; DOA_theta];
         s = s + SNOI;
      otherwise
         error(err_1);
   end;
end;
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%% Thermal Noise %%%%%%%%%%%%%%%%%%%%
if ~isempty(noise)
    for k = 1 : N
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s(k,:) = s(k,:) + noise_data;
    end;
end;
%-------------------------------------------------------%
