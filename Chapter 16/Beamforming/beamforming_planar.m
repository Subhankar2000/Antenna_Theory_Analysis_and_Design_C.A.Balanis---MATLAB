%************************************************************************
%	Beamforming_planar
%************************************************************************
%   It is the Matlab function that simulates beamforming of planar array.
%*************************************************************************
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
[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = planar_data_entry;

% Parameters initialization
FIG           = 'figure(1)';             % Figure to record
w             = zeros(N(1)*N(2),1);

% Generate signals
[dd, X, fm] = planar_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

for i = 1 : length(dd),
    w0 = w;
    [w, err(i)] = LMS(w,Mu,X(:,i),dd(i));
    mse(i)   = sum(abs(err(i))^2);
    w_err(i) = norm(w0 - w);
    if i>1                           
        if 10*log10(mse(i)) < -1000    
            break;
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
disp('=== Iterations Number of the Beamforming ===');              % Displays number of iterations
fprintf('        %d\n',iterationnumber);

fprintf('\n');
disp('===The Exact Weight (Amplitude) of Each Element ===');                 % Displays computed weights in x*y matrix
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d          ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%f   ',W((i-1)*N(2)+j));
    end
    fprintf('\n');
end

fprintf('\n');
disp('===The Normalized Weight (Amplitude) of Each Element ===');                 % Displays computed weights in x*y matrix
fprintf('They are normalized with respect to amplitude of the first element, (x,y)=(1,1).\n');
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d          ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%f   ',W((i-1)*N(2)+j)/W(1));
    end
    fprintf('\n');
end

nbeta=unwrap(beta)*180/pi;
fprintf('\n');
disp('===The Exact Beta (Phase in Degrees) of Each Element ===');                 % Displays computed beta [Ensemble]
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y      ');
for i=1:N(2)
    fprintf('%d             ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%12f  ',nbeta((i-1)*N(2)+j));
    end
    fprintf('\n');
end

fprintf('\n');
disp('===The Normalized Beta (Phase in Degrees) of Each Element ===');                 % Displays computed beta [Ensemble]
fprintf('The phase of the first element, (x,y)=(1,1), is set to 0 degree, and other phases are normalized with repect to it.\n');
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d           ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%10f  ',mod((nbeta((i-1)*N(2)+j)-nbeta(1)),360));
    end
    fprintf('\n');
end

fprintf('\nThe end of the data display. Please have a look at the figures\n');
disp('*************************************************************************');
fprintf('\n');

% Stop recording
if (save_in_file=='y')
    diary off;
end

% Plot results
polar_plot(N,d,w,AF_thresh,E_pattern);
title(['Planar Array Beamforming Pattern (X*Y=', num2str(N(1)),'*',num2str(N(2)),', dx=', num2str(d(1)),'\lambda, dy=', num2str(d(2)), '\lambda, SOI \theta=', num2str(sig(1,2)),'\circ,\phi=',num2str(sig(1,3)), '\circ)']);

AF=planar_plot_pattern(N,d,w,E_pattern,AF_thresh,'planar');
title(['Planar Array Beamforming Pattern (X*Y=', num2str(N(1)),'*',num2str(N(2)),', dx=', num2str(d(1)),'\lambda, dy=', num2str(d(2)), '\lambda, SOI \theta=', num2str(sig(1,2)),'\circ,\phi=',num2str(sig(1,3)), '\circ)']);

figure; plot(W/W(1),'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Excitation of Antenna Element');
title('Normalized Magnitude Distribution');

figure; plot(unwrap(beta')*180/pi,'*-'); grid;
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

% zoom in to see the beamforming pattern in azimuthal plane or elevation
% plane
zoom_in = input('\nDo you want to view the beamforming pattern in azimuthal plane or elevation plane? ([n]/y): ','s');
ERR=1;
while(ERR ~= 0)
   if (isempty(zoom_in) | zoom_in == 'n')
       fprintf('OK. That is all for this simulations.\n\n');
       ERR=0;
   elseif(zoom_in =='y')
       ERR1 = 1;
       while(ERR1 ~= 0)
           fprintf('Plane option: \n\tOption (1): Azimuth\n\tOption (2): Elevation \n');
           plane1 = str2num(input('\n Choice of plane = ','s'));
           if(plane1 == 1)
               fprintf('\n Azimuthal plane is chosen.\n');
               zoom_ele=input('Please enter the value of theta in degrees (between 0 and 180):  ');
               zoom_ele0=zoom_ele;
               if zoom_ele>90    
                   zoom_ele=180-zoom_ele; %simetry of the theta=90.
               end
               figure; 
               plot(linspace(0,360,361),AF(round(zoom_ele*4+1),:));
               grid;
               title(['Beamforming pattern vs \phi at \theta=', num2str(zoom_ele0),'\circ']);
               xlabel('Azimuth (degrees)');
               ylabel('Beamforming pattern (dB)');
           elseif(plane1 == 2)
               fprintf('\n Elevation plane is chosen.\n');
               zoom_azm=input('Please enter the value of phi in degrees (between 0 and 360):  ');
               for i=1:360
                   AF_azm_flip(i)=AF(361-i,round(zoom_azm+1));
               end
               AF_azm=[AF(:,round(zoom_azm+1))',AF_azm_flip];
               figure; 
               plot(linspace(0,180,721),AF_azm);
               grid;
               title(['Beamforming pattern vs \theta at \phi=', num2str(zoom_azm),'\circ']);
               xlabel('Elevation (degrees)');
               ylabel('Beamforming pattern (dB)');
           else
               fprintf('\n Plane number should be either 1 or 2\n');    
           end
               
           zoom_again=input('\n Do you want to try again? ([n]/y): ','s');
           ERR2=1;
           while (ERR2~=0)
               if (isempty(zoom_again) | zoom_again == 'n')
                   fprintf('OK. That is all for this simulations.\n\n')
                   ERR = 0;
                   ERR1 = 0;
                   ERR2 = 0;
               elseif(zoom_again =='y')
                   ERR = 1;
                   ERR1 = 1;
                   ERR2 = 0;
               else
                   zoom_again=input('\n Do you want to try again or not? ([n]/y): ','s');
               end
           end 
       end
   else
      zoom_in = input('\n Do you want to view the beamforming pattern in a specific plane or not? ([n]/y):','s');
   end;
end

% Stop recording
if (save_in_file=='y')
    fprintf('Data saved in file %s.txt\n\n', file_string);
end
