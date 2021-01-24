%************************************************************************
%	Planar.m
%************************************************************************
%   It is the Matlab function that estimate DOAs of planar array.
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

% Set output format data
format short e

% Generate output filename based on current time
cur_time = datevec(datestr(now));

ERR = 1;

save_in_file = input('\nDo you want to save reusults in file? ([n]/y):','s');
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
[N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry;

% Generate signals and estimate the DOAs
DOAs = planar_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
disp('*************************************************************************');
fprintf('\n');
disp('=== Signal Infomation ===');
disp('       Amplitude   Theta(degrees)    Phi(degrees)'); 
[rowofsig,colofsig]=size(sig);
fprintf('SOI     %2.4f        %3d               %3d\n',sig(1,1),sig(1,2),sig(1,3));
for i=2:rowofsig
    fprintf('SNOI%d   %2.4f        %3d               %3d\n',i-1,sig(i,1),sig(i,2),sig(i,3));
end

if(isempty(noise))
    fprintf('\n');
    disp('NO noise in the system.');
else
    fprintf('\n');
    disp('=== Noise Information ===');
    fprintf('Mean:  %f,  Variance:  %f\n', noise(1), noise(2));
end

fprintf('\n');
disp('=== DOA Estimations ===');
disp('        Theta(degrees)  Phi(degrees)');
fprintf('SOI       %3.3f            %3.3f\n', 180/pi*DOAs(1,1),180/pi*DOAs(1,2));
for i=2:rowofsig
    fprintf('SNOI%d     %3.3f            %3.3f\n',i-1,180/pi*DOAs(i,1),180/pi*DOAs(i,2));
end
disp('*************************************************************************');

% Stop recording
if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end
