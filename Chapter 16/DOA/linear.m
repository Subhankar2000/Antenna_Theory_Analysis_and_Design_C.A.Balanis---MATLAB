%*************************************************************************
%	Linear.m
%************************************************************************
%	It is a MATLAB function that estimate the DOA of linear array. 
%************************************************************************
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
[N,d,sig,noise,type,nn,NN,E_pattern] = linear_data_entry;

% Generate signals and estimate DOA.
DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

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
disp('=== DOA Estimations ===');
string = sprintf('        Theta(degrees)'); disp(string);
fprintf('SOI       %3.3f\n', 180/pi*DOAs(1));
for i=2:rowofsig
    fprintf('SNOI%d     %3.3f\n',i-1,180/pi*DOAs(i));
end

disp('*************************************************************************');
% Stop recording

if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end