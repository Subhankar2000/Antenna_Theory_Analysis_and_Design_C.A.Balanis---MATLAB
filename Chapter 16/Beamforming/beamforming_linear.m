%***********************************************************************
%	Beamforming_linear.m
%***********************************************************************
%	It is a MATLAB function that simulates beamforming for linear arrays.
%***********************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

% Set output data format
format short e

% Generate output filename based on current time
cur_time = datevec(datestr(now));

ERR = 1;
save_in_file = input('\n Do you want to save reusults in file? ([n]/y):','s');
while(ERR ~= 0)
   if (isempty(save_in_file) | save_in_file == 'n')
       fprintf('------------Results will not be saved in file.-----------\n\n');
       ERR=0;
   elseif(save_in_file =='y')
       file_string = input('Input the desired output filename: ','s');
       ERR=0;
       fprintf('\n');
   else
       save_in_file=input('Do you want to save the results in file or not? ([n]/y):','s');
   end;
end



% Start timer
tic;

% Start recording
if (save_in_file=='y')
    out_file = sprintf('%s.txt',file_string);
    diary(out_file);
end
    

% User input
[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry;

% Parameters initialization
FIG           = 'figure(1)';             % Figure to record
SKIP_STEP     = 40;                      % Plot every SKIP_STEP iterations
w             = zeros(N,1);              % iteration initialization

% Generate signals
[dd, X, fm] = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

for i = 1 : length(dd)
    w0 = w;
    [w, err(i)] = LMS(w,Mu,X(:,i),dd(i));
    mse(i) = sum(abs(err(i))^2);
    w_err(i) = norm(w0 - w);
    if i>1
        array_factor = linear_AF(N,d,w,sig(1:size(sig,1),2),E_pattern);
        if (abs(w_err(i) - w_err(i-1)) < eps) | (array_factor <= AF_thresh)
            linear_plot_pattern(sig,w,N,d,E_pattern,AF_thresh,'half',4,'-');
            break;
        end;
    end;
    if rem(i,SKIP_STEP) == 0             % Plot every SKIP_STEP iterations
        linear_plot_pattern(sig,w,N,d,E_pattern,AF_thresh,'half',4,'-');
        if i == SKIP_STEP
            FIG_HANDLE = eval(FIG);
        if (isunix)
           pause;
        end;
        end;
    end;
end;
    
% Final weights and betas
W    = abs(w);
beta = angle(w);
iterationnumber=i;

%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
disp('*************************************************************************');
fprintf('\n');
disp('=== Signal Infomation ===');
string = sprintf('       Amplitude   Theta(degrees)'); disp(string);

[rowofsig,colofsig]=size(sig);
fprintf('SOI     %2.4f        %3d\n',sig(1,1),sig(1,2));
for i=2:rowofsig
    fprintf('SNOI%d   %2.4f        %3d\n',i-1,sig(i,1),sig(i,2));
end

if(isempty(noise))
    fprintf('\n');
    disp('NO noise.');
else
    fprintf('\n');
    disp('=== Noise Information ===');
    fprintf('Mean:  %f,  Variance:  %f\n', noise(1), noise(2));
end

fprintf('\n');
disp('=== Iterations Number of Beamforming ===');              % Displays number of iterations
fprintf('        %d\n',iterationnumber);

fprintf('\n');
disp('=== Weights (Amplitude) of Each Element ===');                 % Displays computed weights
disp('   Exaxt Value        Normalized');
for i=1:length(W)
    fprintf('%d   %f           %f\n',i,W(i),W(i)/W(1));
end

fprintf('\n');
disp('=== Beta (Phase in degrees) of Each Element ===');             % Displays computed beta [Ensemble]
disp('      Exaxt Value          Normalized');
nbeta=unwrap(beta)*180/pi;
for i=1:length(beta)
    fprintf('%d   %12f         %12f\n',i,nbeta(i), mod((nbeta(i)-nbeta(1)),360));
end

fprintf('\nThe end.\n');
disp('*************************************************************************');
fprintf('\n');

% Stop recording

if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end


% Plot results
figure; plot(W/W(1),'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Excitation of Antenna Element');
title('Normalized Magnitude Distribution');

figure; plot(unwrap(beta)*180/pi,'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Phase of Antenna Element (degrees)');
title('Phase Distribution');

figure; plot(10*log10(mse)); grid;
xlabel('Iteration Number');
ylabel('MSE in dB');
title('Learning Curve');
    
figure; plot(10*log10(w_err)); grid;
xlabel('Iteration Number');
ylabel('Weight Error in dB');
title('Weight Estimation Error');

% Stop timer
disp('Elapsed Time (min):'); disp(toc/60);

