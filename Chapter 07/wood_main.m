% this program synthesizes a line source or a linear array using the 
% Woodward-Lawson method

% Written by Marios Gkatzianas	
% Arizona State University, September 2002

function []=wood_main(pos);

% load bal.mat;
[x,map]=imread('bal.tiff');

close all;
warning off;

wood_mode=[];
while ~isreal(wood_mode)|isempty(wood_mode)|((wood_mode~=1)&(wood_mode~=2)),
    wood_mode=input( ['Select method of Woodward-Lawson synthesis\n' ...
                      '******************************************\n' ...
                      '1. Line source\n','2. Linear array\n','->'] );
    if ~isreal(wood_mode)|isempty(wood_mode)|(wood_mode~=1&wood_mode~=2),
       hwomod=msgbox('Specifed number must be either 1 or 2','Invalid choice','custom',x,map,'modal');            
    end;               
end;

output_mode=[];
while ~isreal(output_mode)|isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
    output_mode=input([ 'Select output method\n','********************\n', ...
          '1. Screen\n', '2. Output file\n', '->']);
    if ~isreal(output_mode)|isempty(output_mode)|(output_mode~=1&output_mode~=2),
       houtmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
    end;   
end;

if (wood_mode==1),
   str_mod='Space Factor';
else
   str_mod='Array Factor';
end;

outfile=[]; outpath=[];

%%% if (output_mode==2),
%%%	 while isempty(outfile)|(outfile==0), 

%%% don't comment out the next two lines

%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%       if exist('hout'),
%%%           waitfor(hout);
%%%       end;  
%%%       [outfile,outpath]=uiputfile('*.txt','Select output file');
 
%%%       if isa(outfile,'double')|isa(outpath,'double'),     
%%%           delete(gco);                    
%%%           hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                       'custom',x,map,'modal');                 
%%%           set(hout,'interruptible','off','busyaction','queue');
%%%        end;

%%%     end;             
%%% end;

if output_mode==2,
   outfile=input(['Give name of output file to save ' str_mod ' (don''t use quotes)\n'],'s');
   while isempty(outfile),
      outfile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
   end;      
   outpath=pwd;
   outfile=[outpath '\' outfile];   
end;


switch wood_mode
    case 1       
	   len=[];
	   while ~isreal(len)|isempty(len)|len<=0,
         len=input('Give length of line source (in wavelengths)\n');
         if ~isreal(len)|isempty(len)|len<=0,
            hlen=msgbox('Length must be real positive','Invalid choice','custom',x,map,'modal');
         end;  
	   end;

	   M=round(len);
      d=[];
      Nel=[];

    case 2
        Nel=[];
        while ~isreal(Nel)|isempty(Nel)|(Nel<=0)|(mod(real(Nel),1)~=0),
           Nel=input('Specify number of elements for linear array\n');            
           if ~isreal(Nel)|isempty(Nel)|Nel<=0|mod(real(Nel),1)~=0,
              hNel=msgbox('Specified number must be positive integer','Invalid choice', ...
                          'custom',x,map,'modal');
           end;   
        end;
        
	     d=[];
        while ~isreal(d)|isempty(d)|(d<=0),
           d=input('Specify spacing between the elements (in wavelengths)\n');
           if ~isreal(d)|isempty(d)|d<=0,
              hd=msgbox('Specified number must be real positive','Invalid choice','custom',x,map,'modal');
           end;
        end;
        if d>0.5,
           disp(['Element spacing is greater that lambda/2. Insufficient spectral resolution ' ...
                    '(excessive aliasing expexted).']); disp(' ');
        end;
	     len=Nel*d;
        M=round(len);
end;  % switch

sfaf_mode=[];
while ~isreal(sfaf_mode)|isempty(sfaf_mode)|((sfaf_mode~=1)&(sfaf_mode~=2)),
   sfaf_mode=input([str_mod,' definition options\n','********************************\n', ...
         '1. Read ', str_mod,' from a file (angles must be in degrees and ' str_mod ...
         ' in linear scale)\n',...
         '(Note: SF function cannot be explicitly expressed as a function of theta. For example: A rectangular pulse.)\n','\n' ...
         '2. Define ', str_mod,' as a function\n', ...
         '(Note: SF function can be explicitly expressed as a function of theta.)\n', ...
         'Choose one of the above\n->']);
   if ~isreal(sfaf_mode)|isempty(sfaf_mode)|(sfaf_mode~=1&sfaf_mode~=2),
      hsfafmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
   end;   
end;

switch sfaf_mode,
    case 1  % read SF/AF from file. Angles must be in degrees, abs(SF/AF) in 
	    % absolute units   
%        filename=[];
%        while isempty(filename),	 	
%           filename=input(['Specify the name of the file containing ', ...
%                            str_mod,' (don''t use quotes)\n'],'s');
%        end;

        inpath=[]; infile=[];
               
%%%        while isempty(infile)|isa(infile,'double'),
%%%           if exist('hout'),
%%%              waitfor(hout);
%%%           end;   
%%%           [infile,inpath]=uigetfile('*.txt',['Select text file containing ' str_mod]);
%%%           if isa(inpath,'double')|isa(infile,'double'),  
%%%              hout=msgbox('Incorrect file specification','Invalid choice','custom',x,map,'modal');
%%%           end; 
%%%         end;

        infile=input(['Give name of file to read ' str_mod ' (don''t use quotes)\n', ...
              '(Note: The file should have two lines. The first line should contain the values of theta (in degrees),\n', ...
              'and the second line should have the values of SF.\n', ...
              'Read the file ''sfinput.m'' as an example how to create the SF file.\n' ...
              '(''sf1.m'' for rectangular ',str_mod, ','' sf2.m'' for triangular ', str_mod,'))\n'],'s');
        while isempty(infile),
           infile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
        end;      
        inpath=pwd;
        filename=[inpath '\' infile];         
        
        ftid=fopen(filename,'rt');
        if ftid==-1,
           disp('File does not exist or cannot be opened. Program terminated.');
           return;
        else   
           fclose(ftid);
        end;   
           
        SFAF=load(filename);		
%         SFAF=SFAF';
               
        theta=SFAF(1,:);  % degrees
        SFAF=SFAF(2,:);	
	
	     Ntheta=length(theta);	    
	     pstring=[];	          
        
    case 2  

        pstring=[];
        while isempty(pstring),
           pstring=input(['Give ',str_mod,' as a function of x (x=theta).' ...
              ' Use ^ for all powers and don''t use quotes.\n' ...
             '-------------------------------------\n' ... 
             'For example: sin(x)\n' ...
             '-------------------------------------\n'],'s');
           if isempty(pstring),
              hpstr=msgbox('String must be non-empty and contain only a x variable','Invalid choice', ...
                           'custom',x,map,'modal');
           end;              
        end;

        px=inline(vectorize(pstring));

        Ntheta=[];
        while ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
           Ntheta=input(['Specify number of sample points for the angle theta' ...
                         ' for the pattern plots\n', ...
                       '(NOT the sample points for the pattern) [ENTER for default=360]:']);
           Ntheta(isempty(Ntheta))=360;	
           if ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
              hNthe=msgbox('Specified number must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');
           end;   
        end;

	     theta=linspace(0,360,2*Ntheta);
        SFAF=px(theta*pi/180);

	filename=[];

end;   % end of switch

eo_switch=[];
while ~isreal(eo_switch)|isempty(eo_switch)|((eo_switch~=1)&(eo_switch~=2)&(eo_switch~=3)),
    eo_switch=input( ['Select even and/or odd samples\n' ...
                      '******************************************\n' ...
                      '1. Even samples\n','2. Odd samples\n','3. Even & Odd samples\n','->'] );
    if ~isreal(eo_switch)|isempty(eo_switch)|(eo_switch~=1&eo_switch~=2&eo_switch~=3),
       hwomod=msgbox('Specifed number must be either 1,2 or 3','Invalid choice','custom',x,map,'modal');            
    end;               
end;

[theta_even,theta_samp1,b1,z1,I1,SFAF_rec_even]=wood(str_mod,d,len,Nel,sfaf_mode, ...
     filename,pstring,2*M,Ntheta,theta,SFAF);

[theta_odd,theta_samp2,b2,z2,I2,SFAF_rec_odd]=wood(str_mod,d,len,Nel,sfaf_mode, ...
     filename,pstring,2*M+1,Ntheta,theta,SFAF);

U=(abs(SFAF(theta<=180))).^2;  U(isnan(U))=0;
U_even=abs(SFAF_rec_even(theta<=180)).^2; U_even(isnan(U_even))=0;
U_odd=abs(SFAF_rec_odd(theta<=180)).^2; U_odd(isnan(U_odd))=0;

Prad=sum(U.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;
Prad_even=sum(U_even.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;
Prad_odd=sum(U_odd.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;

D=2*max(U)/Prad;
D_even=2*max(U_even)/Prad_even;
D_odd=2*max(U_odd)/Prad_odd;

% Figure 1
% ********
figure(1);

switch eo_switch,
    case 1 % Even samples
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
 	     plot(z1,abs(I1),'k','linewidth',2); grid; hold on; 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['even (' num2str(2*M),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
 	     hs1=plot(z1,abs(I1),'k--','marker','s'); grid; hold on; 
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs1],'linewidth',2,'markersize',6);
        set(hs1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
 	     hl=legend([hs1(1)],['even (' num2str(2*M),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
        
    case 2 % Odd samples
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
	     plot(z2,abs(I2),'linewidth',2,'color','r'); 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['odd (' num2str(2*M+1),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
	     hs2=plot(z2,abs(I2),'r--','marker','d'); grid;
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs2],'linewidth',2,'markersize',6);
        set(hs2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
 	     hl=legend([hs2(1)],['odd (' num2str(2*M+1),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
    case 3 % Both
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
 	     plot(z1,abs(I1),'k','linewidth',2); grid; hold on; 
	     plot(z2,abs(I2),'linewidth',2,'color','r'); 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['even (' num2str(2*M),') samples'],['odd (' num2str(2*M+1),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
 	     hs1=plot(z1,abs(I1),'k--','marker','s'); grid; hold on; 
	     hs2=plot(z2,abs(I2),'r--','marker','d'); grid;
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs1 hs2],'linewidth',2,'markersize',6);
        set(hs1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
        set(hs2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
 	     hl=legend([hs1(1) hs2(1)],['even (' num2str(2*M),') samples'], ...
                   ['odd (' num2str(2*M+1),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
end
grid on;
xlabel('z^{\prime}/\lambda (normalized)','fontname','timesnewroman','fontsize',18, ... 
       'verticalalign','top'); 
ylabel('|I(z^{\prime})|','fontname', 'timesnewroman','fontsize',18,'verticalalign','bottom');

% Figure 2
% ********
figure(2);

switch eo_switch,
    case 1 % Even samples
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_even,abs(SFAF_rec_even),'k','linewidth',2); 
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
hl=legend('Desired',['even (' num2str(2*M),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, even (',num2str(2*M),') samples']);
   case 2
       hl=legend('Desired',['Linear array, even (',num2str(2*M),') samples']);
end
list1=findall(hl,'type','line','color','k');
set(list1,'marker','s','markersize',8,'markerfacecolor','k','markeredgecolor','k');
h1=plot(theta_samp1,abs(b1),'marker','s','linestyle','none');
set(h1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
    case 2 % Odd samples
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_odd,abs(SFAF_rec_odd),'linewidth',2,'color','r');
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
hl=legend('Desired',['odd (' num2str(2*M+1),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, odd (',num2str(2*M+1),') samples']);
   case 2
       hl=legend('Desired',['Linear array, odd (',num2str(2*M+1),') samples']);
end
list2=findall(hl,'type','line','color','r');
set(list2,'marker','d','markersize',8,'markerfacecolor','r','markeredgecolor','r');
h2=plot(theta_samp2,abs(b2),'marker','d','linestyle','none');
set(h2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
    case 3 % Both
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_even,abs(SFAF_rec_even),'k','linewidth',2); 
plot(theta_odd,abs(SFAF_rec_odd),'linewidth',2,'color','r');
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
hl=legend('Desired',['even (' num2str(2*M),') samples'], ...
           ['odd (' num2str(2*M+1),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, even (',num2str(2*M),') samples'],['Line-source, odd (',num2str(2*M+1),') samples']);
   case 2
       hl=legend('Desired',['Linear array, even (',num2str(2*M),') samples'],['Linear array, odd (',num2str(2*M+1),') samples']);
end
list1=findall(hl,'type','line','color','k');
list2=findall(hl,'type','line','color','r');
set(list1,'marker','s','markersize',8,'markerfacecolor','k','markeredgecolor','k');
set(list2,'marker','d','markersize',8,'markerfacecolor','r','markeredgecolor','r');
h1=plot(theta_samp1,abs(b1),'marker','s','linestyle','none');
h2=plot(theta_samp2,abs(b2),'marker','d','linestyle','none');
set(h1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
set(h2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
      
end % end switch

xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','top'); 
ylabel(['|' str_mod '|'],'fontname','timesnewroman','fontsize',18,'verticalalign','bottom');
grid on;


if (output_mode==2),

    fid=fopen(outfile,'wt');
     
    if (wood_mode==1),
       
       fidsaf=fopen('SF_wood.txt','wt');
       fprintf(fid,'%% Woodward-Lawson synthesis method for line source\n\n');
       fprintf(fid,['%% Length of line source = ',num2str(len) ' lambda\n\n']); 
       fprintf(fid,'%% Desired Directivity (dimensionless) = %6.2f\n',D);
       fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Space Factor = %6.2f (even samples)\n',D_even);     
       fprintf(fid,['%%' blanks(34) '>>' blanks(35) '= %6.2f (odd samples)\n\n'],D_odd);         
       fprintf(fid,'%% Desired Directivity (dB) = %6.2f\n',10*log10(D));
       fprintf(fid,'%% Directicity (dB) as computed from synthesized Space Factor = %6.2f (even samples)\n',10*log10(D_even));
       fprintf(fid,['%%' blanks(29) '>>' blanks(29) '= %6.2f (odd samples)\n'],10*log10(D_odd));
                     
       fprintf(fidsaf,'%% Woodward-Lawson synthesis method for line source\n\n');                           
       fprintf(fidsaf,['%% Length of line source = ',num2str(len) ' lambda\n\n']);        
       fprintf(fidsaf,['%%   Theta (deg.)      SF_des      SF_rec (even)' ...
              '      SF_rec (odd)\n']);
       fprintf(fidsaf,'    %8.4f \t     %8.4f \t     %8.4f \t     %8.4f \n',[theta; ...
               SFAF; SFAF_rec_even; SFAF_rec_odd]);
       fclose(fid);	            
    else
       fidsaf=fopen('AF_wood.txt','wt');
                     
       fprintf(fid,'%% Woodward-Lawson synthesis method for linear array\n\n');
       fprintf(fid,['%% Number of array elements = ' num2str(Nel) '\n']);
       fprintf(fid,['%% Spacing = ' num2str(d) ' lambda\n\n']);
       fprintf(fid,'%% Excitation coefficients:\n');              
       fprintf(fid,['%% Real part (' num2str(2*M+1) ' pts) \t Imag part (' num2str(2*M+1) ' pts) \t' ...
                    ' Real part (' num2str(2*M) ' pts) \t Imag part (' num2str(2*M) ' pts)\n']);
       fprintf(fid,'%%   %8.4f \t\t    %8.4f \t\t    %8.4f \t\t   %8.4f\n', ...
               [real(I2); imag(I2); real(I1); imag(I1)]);                    
       fprintf(fid,'\n%% Desired Directivity (dimensionless) = %6.2f\n',D);
       fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Array Factor = %6.2f (even samples)\n',D_even);     
       fprintf(fid,['%%' blanks(34) '>>' blanks(35) '= %6.2f (odd samples)\n\n'],D_odd);         
       fprintf(fid,'%% Desired Directivity (dB) = %6.2f\n',10*log10(D));
       fprintf(fid,'%% Directicity (dB) as computed from synthesized Array Factor = %6.2f (even samples)\n',10*log10(D_even));
       fprintf(fid,['%%' blanks(29) '>>' blanks(29) '= %6.2f (odd samples)\n'],10*log10(D_odd));
              
       fprintf(fidsaf,'%% Woodward-Lawson synthesis method for linear array\n\n');
       fprintf(fidsaf,['%% Number of array elements = ' num2str(Nel) '\n']);
       fprintf(fidsaf,['%% Spacing = ' num2str(d) ' lambda\n\n']);       
       fprintf(fidsaf,['\n%%   Theta (deg.)      AF_des      AF_rec (even)' ...
              '      AF_rec (odd)\n']);
       fprintf(fidsaf,'    %8.4f \t     %8.4f \t     %8.4f \t     %8.4f \n',[theta; ...
               SFAF; SFAF_rec_even; SFAF_rec_odd]);
       fclose(fid);     
    end;
    
end;

disp(' ');
disp(['Desired Directivity (dimensionless) = ' num2str(D,'%6.2f')]);   
disp(['Directivity (dimensionless) as computed from synthesized ' str_mod ' = ' num2str(D_even,'%6.2f') ...
      ' (even samples)']);
disp([blanks(34) '>>' blanks(34) '= ' num2str(D_odd,'%6.2f') ' (odd samples)']);
disp(' ');
disp(['Desired Directivity (dB) = ' num2str(10*log10(D),'%6.2f')]);   
disp(['Directivity (dB) as computed from synthesized ' str_mod ' = ' num2str(10*log10(D_even),'%6.2f') ...
      ' (even samples)']);
disp([blanks(29) '>>' blanks(28) '= ' num2str(10*log10(D_odd),'%6.2f') ' (odd samples)']);
      

if output_mode==2,
   disp(' ');
   disp(['Output data saved in file ' '''' outfile '''']);
   if wood_mode==1,
      disp([str_mod ' saved in file ''' pwd '\SF_wood.txt''']);
   else
      disp([str_mod ' saved in file ''' pwd '\AF_wood.txt''']);
   end;   
end;   

disp(' ');
disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;

fclose('all');



