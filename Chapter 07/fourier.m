% this program synthesizes a line source or a linear array using 
% the Fourier method

% Written by Marios Gkatzianas  
% Arizona State University, September 2002

close all;
clear all;

% load bal.mat;
[x,map]=imread('bal.tiff');

fourier_mode=[];
while isempty(fourier_mode)|((fourier_mode~=1)&(fourier_mode~=2)),
   fourier_mode=input(['Select method of Fourier synthesis\n', ...
                       '**********************************\n', ...
                       '1. Line source\n','2. Linear array\n','->' ]);  
                 
   if isempty(fourier_mode)|((fourier_mode~=1)&(fourier_mode~=2)),
      hsof=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal'); 
   end;                                  
end;

output_mode=[];
while isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
   output_mode=input(['Select output method\n','********************\n', ...
         '1. Screen\n', '2. Output file\n','->' ]);
   
   if isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
      hsouf=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal'); 
   end;   
end;

if fourier_mode==1,
   str_mod='Space Factor';
else
   str_mod='Array Factor';
end;   

outfile=[]; 
%%% if (output_mode==2),
%%%	  while isempty(outfile)|(outfile==0), 

%%% don't comment out the next two lines

%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%           if exist('hout'),
%%%               waitfor(hout);
%%%           end;  
%%%           [outfile,outpath]=uiputfile('*.txt','Save As');
 
%%%           if isa(outfile,'double')|isa(outpath,'double'),     
%%%               delete(gco);                    
%%%               hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                           'custom',x,map,'modal');                 
%%%               set(hout,'interruptible','off','busyaction','queue');
%%%           end;

%%%       end;
             
%%% end;

if output_mode==2,
   while isempty(outfile),
      outfile=input(['Give name of output file to save ' str_mod ' (don''t use quotes)\n'],'s');
   end;
   outpath=pwd;
   outfile=[outpath '\' outfile];
end;

switch fourier_mode,
    case 1
       str_mode='SF';
       len=[];
       while isempty(len)|len<=0,
          len=input('Give length of line source (in wavelengths)\n');   
          
          if isempty(len)|(len<=0),
             hsoff=msgbox('Length must be non-empty and real positive','Invalid choice','custom',x,map,'modal'); 
          end;                    
       end;
                                   
    case 2
       str_mode='AF';
       Nel=[];
       while isempty(Nel)|Nel<=0|mod(Nel,2)~=1,
          Nel=input(['Give number of elements for the array. Must be an odd number since there ' ,...
                     'is a non-zero dc value.\n']);
          
          if isempty(Nel)|Nel<=0|mod(Nel,2)~=1,
             hsnf=msgbox('Specified number must be positive odd integer','Invalid choice','custom',x,map,'modal'); 
          end;                               
       end;

       d=[];

       disp('It is assumed that d=0.5*lambda and beta=0'); disp(' ');
       d=0.5;

end;  % switch


sfaf_mode=[];
while isempty(sfaf_mode)|((sfaf_mode~=1)&(sfaf_mode~=2)),
   sfaf_mode=input([str_mode,' definition options\n', ...
         '*********************\n','1. Read ', str_mode,' from file (angles must be in degrees and ', ...
         str_mode,' in linear scale)\n', ...
         '(Note: SF function cannot be explicitly expressed as a function of theta. For example: A rectangular pulse.)\n','\n' ...
      '2. Specify ',str_mode,' as a function\n', ...
      '(Note: SF function can be explicitly expressed as a function of theta.)\n', ...
      'Choose one of the above\n','->' ]);

   if isempty(sfaf_mode)|((sfaf_mode~=1)&(sfaf_mode~=2)),
      hsff=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal'); 
   end;  
end;

%%% if sfaf_mode==1,

%%%   inpath=[];
%%%   infile=[];
               
%%%   while isempty(infile)|isa(infile,'double'),
%%%       if exist('hout'),
%%%           waitfor(hout);
%%%       end;   
%%%       [infile,inpath]=uigetfile('*.txt',['Select text file to read ' str_mode ' from']);
%%%       if isa(inpath,'double')|isa(infile,'double'),  
%%%           hout=msgbox('Incorrect file specification','Invalid choice','custom',x,map,'modal');
%%%       end; 
%%%   end;
    
%%%   inname=[inpath infile];
%%% end;

if sfaf_mode==1,
   inpath=[]; infile=[];
   
   while isempty(infile),
      infile=input(['Give name of file to read ' str_mod ' (don''t use quotes)\n', ...
              '(Note: The file should have two lines. The first line should contain the values of theta (in degrees)\n',.... 
              'and the second line should have the values of SF).\n', ...
              'Read the file ''sfinput.m'' as an example how to create the SF file.\n' ...
              '(''sf1.m'' for rectangular ',str_mod, ','' sf2.m'' for triangular ', str_mod,'))\n'],'s');
   end;   
   inpath=pwd;
   infile=[inpath '\' infile];
   ftid=fopen(infile);
   if ftid==-1,
      disp('File does not exist or cannot be accessed. Program terminated.');
   end;   
   fclose(ftid);
end;

   
switch sfaf_mode,

    case 1    % read SF from file. It is assumed that angles are in degrees
       SFAF=load(infile);
       theta=SFAF(1,:)*pi/180;  % convert to radians
       Ntheta=length(theta);
       SFAF=abs(SFAF(2,:));    % amplitude is taken to protect the user
       
       st=chap(SFAF,str_mod);
%%%       if st==1,
%%%          return;
%%%       end;   
       
    case 2    % define SF through inline function
   
       Ntheta=[];
       while isempty(Ntheta)|Ntheta<=0|mod(Ntheta,1)~=0,
          Ntheta=input(['Specify number of samples for the angle theta for the pattern plots ' ...
                       '[ENTER for default=360]:']);
          Ntheta(isempty(Ntheta))=360;
       end;
       
       theta=linspace(0,pi,Ntheta);

       fstring=[];
       while isempty(fstring),
          fstring=input(['Give ',str_mode,' as a function of x (x=theta).' ...
              ' Use ^ for all powers and don''t use quotes.\n' ...
             '-------------------------------------\n' ... 
             'For example: sin(x)\n' ...
             '-------------------------------------\n'],'s');
       end;
       
       ffun=inline(vectorize(fstring));
       SFAF=abs(ffun(theta));             
       
       st=chap(SFAF,str_mod);
%%%       if st==1,
%%%          return;
%%%       end;   
       
end; % switch      


if str_mode=='AF',

   psi=pi*cos(theta);
   dth=theta(2)-theta(1);

   if mod(Nel,2)~=0,    % odd samples

     M=(Nel-1)/2;
     m=-M:M;

     for ind=1:length(m),
        a(ind)=1/2*sum(SFAF.*exp(-j*m(ind)*pi*cos(theta)).*sin(theta)*dth);
     end;
     
     SFAF_rec=zeros(size(theta));
     
     for ind=1:length(m),
        SFAF_rec=SFAF_rec+a(ind)*exp(j*m(ind)*psi);
     end;               
     
   end;  % odd samples
   
   if mod(Nel,2)==0,   % even samples
   
      M=Nel/2;
      m=[-M:-1 1:M];
   
      for ind=1:length(m),      
         if m(ind)>0,
            a(ind)=1/2*sum(SFAF.*exp(-j*(2*m(ind)+1)/2*psi).*sin(theta)*dth);
         else   
            a(ind)=1/2*sum(SFAF.*exp(-j*(2*m(ind)-1)/2*psi).*sin(theta)*dth);
         end;            
      end;
      
      SFAF_rec=zeros(size(theta));
      for ind=1:length(m),
         SFAF_rec=SFAF_rec+a(ind)*exp(j*m(ind)*psi);
      end; 
   end; % even samples

end;  % arrays


if str_mode=='SF',   % line source

   xi=2*pi*cos(theta);
   zs=linspace(-len/2,len/2,Ntheta);
   dzs=zs(2)-zs(1);
   if dzs>1/(4*pi),
      disp('Insufficient number of sample points. Program terminated');
      return;
   end;

   dth=theta(2)-theta(1);

%  if abs(SFAF(1))<=1e-3,   % compact support                   
      
   for n=1:length(zs),
       Is(n)=sum(SFAF.*exp(-j*2*pi*zs(n)*cos(theta)).*sin(theta)*dth);
   end;
      
   for n=1:Ntheta, 
       SFAF_rec(n)=sum(Is.*exp(j*2*pi*zs*cos(theta(n)))*dzs);         
   end;
      
%   end;   
      
   if abs(SFAF(1))>=1e-3, % expand xie domain with exponential decay in invisible region
      
      invi=2.1;
      xie=linspace(-2*pi,2*pi,Ntheta);      
      dxie=xie(2)-xie(1);
      xiel=-2*pi*invi:dxie:(-2*pi-dxie);
      xier=(2*pi+dxie):dxie:2*pi*invi;
      expl=SFAF(end)*exp(-5/((invi-1)*pi)*abs(xiel+2*pi));
      expr=SFAF(1)*exp(-5/((invi-1)*pi)*abs(xier-2*pi));
      
      xie=[xiel xie xier];
      SFAF=[expl fliplr(SFAF) expr];     
      
      zs2=linspace(-len/2,len/2,length(SFAF));
      dzs=zs2(2)-zs2(1);
      
      if dzs>1/(4*pi),
         disp('Insufficient number of sample points. Program terminated');
         return;
      end;           
            
      for n=1:length(zs2),
         Is2(n)=1/(2*pi)*sum(SFAF.*exp(-j*xie*zs2(n))*dxie);
      end;
      
      for n=1:length(theta), 
         SFAF_rec2(n)=sum(Is2.*exp(j*2*pi*zs2*cos(theta(n)))*dzs);
      end;                               
                           
      
   end;  % if
             
end;  % line source


if exist('SFAF_rec2'),
   U=abs(fliplr(SFAF(length(expl)+1:length(expl)+length(theta)))).^2;   
else
   U=abs(SFAF).^2;
end;   

D=2*max(U(theta<pi))/sum(U(theta<pi).*sin(theta(theta<pi))*dth);
U_rec=abs(SFAF_rec).^2;
D_rec=2*max(U_rec(theta<pi))/sum(U_rec(theta<pi).*sin(theta(theta<pi))*dth);

if exist('SFAF_rec2'),
   U_rec2=abs(SFAF_rec2).^2;
   D_rec2=2*max(U_rec2(theta<pi))/sum(U_rec2(theta<pi).*sin(theta(theta<pi))*dth);
end;   

disp(' ');
disp(['Desired directivity (dimensionless) based on ' str_mod ' = ' num2str(D,'%6.2f')]);
disp(['Desired directivity (dB) based on ' str_mod ' = ' num2str(10*log10(D),'%6.2f')]);
disp(' ');
disp(['Directivity (dimensionless) of reconstructed ' str_mod ' (visible region only) = ' num2str(D_rec,'%6.2f')]);
disp(['Directivity (dB) of reconstructed ' str_mod ' (visible region only) = ' num2str(10*log10(D_rec),'%6.2f')]);

if exist('SFAF_rec2'),
   disp(' ');
   disp(['Directivity (dimensionless) of reconstructed ' str_mod ' (invisible region included) = ' num2str(D_rec2,'%6.2f')]);
   disp(['Directivity (dB) of reconstructed ' str_mod ' (invisible region included) = ' num2str(10*log10(D_rec2),'%6.2f')]);
end;

[pre,HPBW_app,tm]=find_hpbw(theta,U);
[pre1,HPBW_rec,tm1]=find_hpbw(theta,U_rec);

if exist('SFAF_rec2'),
   [pre2,HPBW_rec2,tm2]=find_hpbw(theta,U_rec2);
end;   

disp(' ');
disp(['Maximum radiation of desired ' str_mod ' at ' num2str(tm*180/pi) ' deg.']);
if pre<=1e-2,
   disp(['HPBW of desired ' str_mod ' = ' num2str(HPBW_app,'%6.2f') ' deg.']);
else
   disp(['HPBW of desired ' str_mod ' could not be computed with sufficient accuracy.' ...
         ' You must determine it visually.']);     
end;   

disp(' ');
disp(['Maximum radiation of reconstructed ' str_mod ' (visible region only) at ' num2str(tm1*180/pi) ' deg.']);
if pre1<=1e-2,
   disp(['HPBW of reconstructed ' str_mod ' (visible region only) = ' num2str(HPBW_rec,'%6.2f') ' deg.']);
else   
   disp(['HPBW of reconstructed ' str_mod ' (visible region only) could not be computed with sufficient accuracy.']);
   disp('You must determine it visually.');                
end;

if exist('SFAF_rec2'),
   disp(' ');
   disp(['Maximum radiation of reconstructed ' str_mod ' (invisible region included) at ' num2str(tm2*180/pi) ' deg.']);
   if pre2<=1e-2,
      disp(['HPBW of reconstructed ' str_mod ' (invisible region included) = ' num2str(HPBW_rec2,'%6.2f') ' deg.']);  
   else
      disp(['HPBW of reconstructed ' str_mod ' (invisible region included) could not be determined analytically.']); 
      disp('You must determine it visually.');
   end;   
end;

if (output_mode==2),   
   fid=fopen(outfile,'wt');
   fmid=fopen([str_mode '_fourier.txt'],'wt');   
   fprintf(fid,'%% Fourier synthesis method\n\n'); 
   fprintf(fmid,'%% Fourier synthesis method\n\n');
   
   if str_mode=='SF',   
      fprintf(fid,'%% Length of line source = %6.2f lambda\n\n',len); 
      fprintf(fmid,'%% Length of line source = %6.2f lambda\n\n',len);
   else
      fprintf(fid,'%% Number of elements= %3d\n',Nel);
      fprintf(fmid,'%% Number of elements= %3d\n',Nel);
      fprintf(fid,'%% Spacing d=0.5 lambda\n\n');
      fprintf(fmid,'%% Spacing d=0.5 lambda\n\n');
   end;
   
   
   if str_mode=='SF'&length(SFAF)==length(theta),
      fprintf(fmid,'%% theta (deg.) \t SF_des \t SF_rec\n');
   elseif str_mode=='SF'&length(SFAF)~=length(theta),
      fprintf(fmid,'%% theta (deg.) \t SF_des \t SF_rec (vis only) \t SF_rec (invis)\n');
   else   % array factor
      fprintf(fid,'%% Excitation coefficients:\n');
      fprintf(fid,'%% Real part \t Imag part\n');
      fprintf(fid,'%10.4f %10.4f\n',[real(a); imag(a)]);
      fprintf(fmid,'\n%% theta (deg.) \t AF_des \t AF_rec\n');    
   end;   
   
   if length(SFAF)~=length(theta),
      fprintf(fmid,'%8.4f \t %8.4f \t %8.4f \t %8.4f\n',[theta*180/pi; ...
            abs(SFAF(length(expl)+1:length(expl)+length(theta))); abs(SFAF_rec); abs(SFAF_rec2)]);              
   else         
      fprintf(fmid,'%8.4f \t %8.4f \t %8.4f\n',[theta*180/pi; abs(SFAF); abs(SFAF_rec)]);                 
   end;   
   
   fprintf(fid,['\n%% Desired directivity (dimensionless) based on ' str_mod ' = %6.2f\n'],D);
   fprintf(fid,['%% Desired directivity (dB) based in ' str_mod ' = %6.2f\n\n'],10*log10(D));
   fprintf(fid,['%% Directivity (dimensionless) of reconstructed ' str_mod ' (visible region only) = %6.2f\n'],D_rec);
   fprintf(fid,['%% Directivity (dB) of reconstructed ' str_mod ' (visible region only) = %6.2f\n'],10*log10(D_rec));
   
   if exist('SFAF_rec2'),
      fprintf(fid,['\n%% Directivity (dimensionless) of reconstructed ' str_mod ' (invisible region included) = %6.2f\n'], D_rec2);
      fprintf(fid,['%% Directivity (dB) of reconstructed ' str_mod ' (invisible region included) = %6.2f\n'], 10*log10(D_rec2));                     
   end;         
    
   fprintf(fid,['\n%% Maximum radiation on desired ' str_mod ' at %6.2f deg.\n'],tm*180/pi);
   if pre<=1e-2,
      fprintf(fid,['%% HPBW of desired ' str_mod ' = %6.2f deg.\n'],HPBW_app);
   else
      fprintf(fid,['%% HPBW of desired ' str_mod ' could not be computed with sufficient accuracy.\n']);
      fprintf(fid,'%% You must determine it visually.\n');
   end;   
   
   fprintf(fid,['\n%% Maximum radiation of reconstructed ' str_mod ' (visible region only) at %6.2f deg.\n'],tm1*180/pi);
   if pre1<=1e-2,
      fprintf(fid,['%% HPBW of reconstructed ' str_mod ' (visible region only) = %6.2f deg.\n'],HPBW_rec);
   else      
      fprintf(fid,['%% HPBW of reconstructed ' str_mod ' (visible region only) could not be computed with sufficient accuracy.\n']);
      fprintf(fid,'%% You must determine it visually.\n');      
   end;   
      
   if exist('SFAF_rec2'),
      fprintf(fid,['\n%% Maximum radiation of reconstructed ' str_mod ' (invisible region inlcuded) at %6.2f deg.\n'],tm2*180/pi);      
      if pre2<=1e-2,
         fprintf(fid,['%% HPBW of reconstructed ' str_mod ' (invisible region included) = %6.2f deg.\n'],HPBW_rec2);
      else         
         fprintf(fid,['%% HPBW of reconstructed ' str_mod ' (invisible region included) could not be computed with sufficient accuracy.\n']);
         fprintf(fid,'%% You must determine it visually.\n');         
      end;   
   end;
    
end;     


% Figure 1
% ********
if exist('SFAF_rec2'),  % invisible region    
   plot(zs,abs(Is),'linewidth',2,'color','r'); hold on; plot(zs2,abs(Is2),'linewidth',2,'color','g');
   set(gca,'fontname','timesnewroman','fontsize',14);
   legend('Visible region only','Invisible region included'); grid on;
elseif str_mode=='SF',
   plot(zs,abs(Is),'linewidth',2);
   set(gca,'fontname','timesnewroman','fontsize',14);
end;  

if str_mode=='SF',
  xlabel('z^{\prime}/\lambda','fontname','timesnewroman','fontsize',18);
  ylabel('Amplitude distribution','fontname','timesnewroman','fontsize',18);
  title(['Current amplitude using the F-T method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman','fontsize',18);
 grid on;  
end;

if str_mode=='AF',
   za=(-M:M)*d;
   subplot(2,1,1); hst1=plot(za,abs(a),'--','marker','d','linewidth',2);
   set(gca,'fontname','timesnewroman','fontsize',14);
   ylabel('Amplitude distribution ','fontname','timesnewroman','fontsize',18);
   title(['Amplitude and phase using the F-T method for the array',' (N = ',num2str(Nel),')'],'fontname','timesnewroman', ...
       'fontsize',18);   grid on;  
   subplot(2,1,2); hst2=plot(za,angle(a)*180/pi,'--','marker','d','linewidth',2); grid on;
   set(gca,'fontname','timesnewroman','fontsize',14);
   set([hst1(:) hst2(:)],'linewidth',2);
   xlabel('z^{\prime}/\lambda','fontname','timesnewroman','fontsize',18);
   ylabel('Phase distribution','fontname','timesnewroman','fontsize',18);
   set(hst1,{'markerfacecolor','markeredgecolor','markersize'}, ...
           {'b' 'b' 8});
   set(hst2,{'markerfacecolor','markeredgecolor','markersize'}, ...
           {'b' 'b' 8});      
end;

% Figure 2
% ********
figure(2);

if length(SFAF)~=length(theta),
    plot(theta*180/pi,abs(fliplr(SFAF(length(expl)+1:length(expl)+length(theta)))),'color','b','linewidth',2);
    hold on;
    plot(theta*180/pi,abs(SFAF_rec),'linewidth',2,'color','r'); 
    plot(theta*180/pi,abs(SFAF_rec2),'linewidth',2,'color','g');
    set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
else
    plot(theta*180/pi,abs(SFAF),'color','b','linewidth',2); hold on;
    plot(theta*180/pi,abs(SFAF_rec),'color','r','linewidth',2);  
    set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
end;

if str_mode=='SF',
    legend('Desired',['Line-source (','L = ',num2str(len),'\lambda)']); grid on;
else
    legend('Desired',['Linear array (N = ',num2str(Nel),')']); grid on;
end


xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18);
ylabel(['|' str_mod '|'],'fontname','timesnewroman','fontsize',18,'verticalalign','bottom');
title(['Synthesized ' str_mod ' using the F-T method'], ...
       'fontname','timesnewroman','fontsize',18);
if exist('fid'),
   fclose(fid);
end;

if exist('SFAF_rec2'),
   mse1=mean(abs(SFAF(length(expl)+1:length(expl)+length(theta))-SFAF_rec).^2);
   mse2=mean(abs(SFAF(length(expl)+1:length(expl)+length(theta))-SFAF_rec2).^2);   
   disp(' ');
   disp(['Mean square error of reconstruction (visible region only) = ',num2str(mse1,'%8.4f')]);
   disp(['Mean square error of reconstruction (invisible region included) = ',num2str(mse2,'%8.4f')]);   
end;   

disp(' ');
if output_mode==2,
   disp(['Output data saved in file ''',outfile,'''']);
   disp([str_mod ' saved in file ''' pwd '\' str_mode '_fourier.txt''']);
end;   
disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;
fclose('all');

