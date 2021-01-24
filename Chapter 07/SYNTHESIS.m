%*************************************************************************************************
% SYNTHESIS.m
%*************************************************************************
% This is a MATLAB based program that:
%
% Implements the various synthesis methods presented in Chapter 7 of Antenna Theory and Design. 
% 
% Specifically, the:
%
% A. Schelkunoff
% B. Fourier
% C. Woodward-Lawson
% D. Taylor (both Tschebyscheff-Error and One-Parameter) 
% 
% methods are used to synthesize line sources and/or linear arrays (where applicable). 
%
% The user is guided in entering the correct input data through interactive questions 
% with built-in error checking. 
%
% The output data:
% A. Are plotted
% B. Can be saved both in MATLAB's binary .mat format as well as user-specified text files. 
% .........................................................................
% Written by Marios Gkatzianas
% Arizona State University, September 2002
% .........................................................................
% Revised by Bo Yang
% Arizona State University, September 2004
% ........................................................................

close all;
clear all;

figure(1)
%%% disp('Please read the following carefully:');
%%% disp('The MATLAB prompt is in pause mode right now and an empty MATLAB figure is shown.');
%%% disp('Resize the size of the figure with the mouse (DON''T PRESS ANY KEY) until you are');
%%% disp('satisfied and THEN PRESS any key. The figure size you have selected will be used');
%%% disp('throughout this MATLAB program for all subsequent figures.');
%%% pause;

pos=get(gcf,'position');
close all;

disp(' ');
disp('Figure position acquired.'); disp(' ');
load bal.mat;

form=5;
while (form~=1)&(form~=2)&(form~=3)&(form~=4),
    form=input(['Select format for output data\n','*****************************\n', ...
         '1. short (5 decimal digits)\n','2. long (15 decimal digits)\n','3. short e\n','4. long e\n', ...
         'Select desired format (ENTER for default=short):']);
    form(isempty(form))=1; 
    
    if (form~=1)&(form~=2)&(form~=3)&(form~=4),
       hform=msgbox('Specified number must be either 1,2,3 or 4','Invalid choice','custom',x,map,'modal'); 
    end;                
end;

switch form
   case 1,
      format short;
   case 2
      format long;
   case 3
      format short e;
   otherwise
      format long e;
end;


method=[];
while isempty(method)|((method~=1)&(method~=2)&(method~=3)&(method~=4)),
   method=input(['Choose one of the following synthesis methods\n', ...
          '*********************************************\n' ...
          '1. Schelkunoff method\n', '2. Fourier transform method\n', ...
          '3. Woodward-Lawson method\n', '4. Taylor method\n','->']);
        
    if isempty(method)|((method~=1)&(method~=2)&(method~=3)&(method~=4)),
        hm=msgbox('Specified number must be either 1,2,3 or 4','Invalid choice','custom',x,map,'modal'); 
    end;   
end;

switch method
   case 1,    % Schelkunoff

        output_mode=[];
        while isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
            output_mode=input(['Select output method\n', ...
                     '********************\n', '1. Screen\n', ...
                     '2. Output file\n', '->']);
            if isempty(output_mode)|(output_mode~=1&output_mode~=2);
                hf=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');  
            end;      
        end;		
	
        outfile=[]; outpath=[];
   
%%%	     if (output_mode==2),
%%%	         while isempty(outfile)|(outfile==0), 

% don't comment out the next two lines
%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%                 if exist('hout'),
%%%                    waitfor(hout);
%%%                 end;  
%%%                 [outfile,outpath]=uiputfile('*.txt','Save As');
 
%%%                 if isa(outfile,'double')|isa(outpath,'double'),     
%%%                     delete(gco);                    
%%%                     hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                                 'custom',x,map,'modal');                 
%%%                     set(hout,'interruptible','off','busyaction','queue');
%%%                 end;

%%%             end;
             
%%% 	     end;

        if output_mode==2,
             outfile=input(['Give name of output file to save Array Factor (don''t use quotes)\n'],'s');
             while isempty(outfile),
                outfile=input(['Filename must be non-empty. Re-enter name of output file ' ...
                              '(don''t use quotes).\n'],'s');
             end;     
        end;
   
        outfile=[pwd '\' outfile];
        schel(output_mode,outfile,pos);

   case 2     % Fourier
      
        fourier;

   case 3,    % Woodward

	wood_main(pos);

   case 4,    % Taylor

	taylor_main(pos);

end;  % switch