% this program synthesizes a line source using the Taylor method 
% variations of Chebyshev error and one parameter

% Written by Marios Gkatzianas  
% Arizona State University, September 2002

function []=taylor_main(pos);

% load bal.mat; 
[x,map]=imread('bal.tiff');

close all;

taylor_mode=[];
while ~isreal(taylor_mode)|isempty(taylor_mode)|((taylor_mode~=1)&(taylor_mode~=2)),
   taylor_mode=input( ['Select method of Taylor line source synthesis\n' ...
                     '*********************************************\n' ...
                     '1. Tschebyscheff error\n','2. One-Parameter\n','->'] );
   if ~isreal(taylor_mode)|isempty(taylor_mode)|(taylor_mode~=1&taylor_mode~=2),           
      htaymod=msgbox('Specified number must be either 1 or 2' ,'Invalid choice','custom',x,map,'modal');     
   end;   
end;

output_mode=[];
while ~isreal(output_mode)|isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
    output_mode=input(['Select output method\n','********************\n', ...
          '1. Screen\n','2. Output file\n','->']);
    if ~isreal(output_mode)|isempty(output_mode)|(output_mode~=1&output_mode~=2),
       houtmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
    end;   
end;


outfile=[]; outpath=[];

%%% if (output_mode==2),
%%%	 while isempty(outfile)|(outfile==0), 

% don't comment in the next two lines for GUI style

%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%       if exist('hout'),
%%%           waitfor(hout);
%%%       end;  
%%%       [outfile,outpath]=uiputfile('*.txt','Save As');
 
%%%       if isa(outfile,'double')|isa(outpath,'double'),     
%%%           delete(gco);                    
%%%           hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                       'custom',x,map,'modal');                 
%%%           set(hout,'interruptible','off','busyaction','queue');
%%%        end;

%%%     end;
             
%%% end;
         
%%% outfile=[outpath outfile];

if output_mode==2,
   outfile=input('Give name of output file to save Space Factor (don''t use quotes)\n','s');
   while isempty(outfile),
      outfile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
   end;      
   outpath=pwd;
   outfile=[outpath '\' outfile];   
end;

len=[];
while ~isreal(len)|isempty(len)|len<=0,
   len=input('Specify length of line source (in wavelengths)\n');
   if ~isreal(len)|isempty(len)|len<=0,
      hlen=msgbox('Length must be non-empty and real positive','Invalid choice','custom',x,map,'modal');
   end;   
end;

side_lobe=[];
while ~isreal(side_lobe)|isempty(side_lobe)|side_lobe>=0,
   if taylor_mode==1,
      side_lobe=input('Specify desired constant side lobe level (in -dB)\n');
   else
      side_lobe=input('Specify desired maximum side lobe level (in -dB)\n');
   end;   
   if ~isreal(side_lobe)|isempty(side_lobe)|side_lobe>=0,
      hslobe=msgbox('Side lobe level must be real negative','Invalid choice','custom',x,map,'modal');
   end;   
end;

if (taylor_mode==1),
    nbar=[];
    while ~isreal(nbar)|isempty(nbar)|nbar<=0|mod(nbar,1)~=0|real(nbar)>len, 
       nbar=input('Specify number of side lobes of the same constant level (nbar)\n');
       if ~isreal(nbar)|isempty(nbar)|nbar<=0|mod(real(nbar),1)~=0|real(nbar)>len,          
          hnbar=msgbox(strvcat('nbar must be a positive integer less than the' , ... 
                    'length of the line source'),'Invalid choice','custom',x,map,'modal');
       end;   
    end;
end;

Ntheta=[];
while ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0, 
   Ntheta=input(['Specify number of samples for the angle theta for the pattern plots\n' ...
                 '[ENTER for default=360]:']);
   Ntheta(isempty(Ntheta))=360;
   if ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,        
      hNthe=msgbox('Specified number must be positive integer','Invalid choice','custom',x,map,'modal');
   end;   
end;

theta=linspace(0,pi,Ntheta);
u=pi*len*cos(theta);

Nz=[];
while ~isreal(Nz)|isempty(Nz)|Nz<=0|mod(real(Nz),1)~=0,
   Nz=input(['Specify number of sample points for the line source for plotting purposes\n' ...
            '[ENTER for default=200]:']);
   Nz(isempty(Nz))=200;
   if ~isreal(Nz)|isempty(Nz)|Nz<=0|mod(real(Nz),1)~=0,   
      hNz=msgbox('Specified number must be positive integer','Invalid choice','custom',x,map,'modal');
   end;   
end;

thres=[];
while ~isreal(thres)|isempty(thres)|thres>=0,
   thres=input(['Specify the lowest threshold for dB scale (in -dB) for plotting purposes ' ...
               '[ENTER for default=-120]:']);
   thres(isempty(thres))=-120;
   if ~isreal(thres)|isempty(thres)|thres>=0,
      hthre=msgbox('Specified number must be real negative','Invalid choice','custom',x,map,'modal');
   end;   
end;

R0=10^(-side_lobe/20);

switch taylor_mode 
    case 1   % Chebyshev error

	A=1/pi*acosh(R0);
	sig=nbar/sqrt(A^2+(nbar-0.5)^2);

	n=1:nbar-1;
	costh1=sig/len*sqrt(A^2+(n-0.5).^2);
	un=pi*len*costh1(abs(costh1)<=1);
	costh1=[-costh1 costh1];

	ng=nbar:floor(len);
	costh2=ng/len;
	costh2=[-costh2 costh2];

	costh=[costh1 costh2];
	costh(abs(costh)>1)=[];

	theta_nul=sort(acos(costh))*180/pi;  % side lobes with constant level

	SF_coef=ones(size(un));

	for p=1:nbar-1,
	    SF_coef(p)=(fact(nbar-1))^2/(fact(nbar-1+p)*fact(nbar-1-p))* ...
                        prod(1-(p*pi./un).^2);
	end

	zs=linspace(-len/2,len/2,Nz);
	Is=ones(size(zs));

	for p=1:nbar-1,
	   Is=Is+2*SF_coef(p)*cos(2*pi*p*zs/len);
	end;
	Is=Is/len;

SF_rec=sinc(u/pi);   

   
	for p=1:nbar-1,
	   SF_rec=SF_rec.*(1-(u/un(p)).^2)./(1-(u/(p*pi)).^2);
	end;

     case 2   % one parameter

	     lobe_look=[-10:-5:-40];  % lookup table
	     B=[j*0.4597 0.3558 0.7386 1.0229 1.2761 1.5136 1.7415];
        
        B=B(find(lobe_look==side_lobe));

        if isempty(B),   % find B numerically           
	        f=inline(['sinh(pi*x)/(pi*x)-', num2str(R0), '/4.603']); 
	        x=fzero(f,1,1e-3);
	        B=x;
        end;
  
	     zs=linspace(-len/2,len/2,Nz);
	     Is=besselj(0,j*pi*B*sqrt(1-(2*zs/len).^2));

	     SF_rec=zeros(size(theta));

 	     SF_rec(abs(u)<pi*B)=len*sinh(sqrt((pi*B)^2-u(abs(u)<pi*B).^2))./ ...
                            sqrt((pi*B)^2-u(abs(u)<pi*B).^2);

        SF_rec(abs(u)>pi*B)=len*sin(sqrt(u(abs(u)>pi*B).^2-(pi*B)^2))./ ...
			    sqrt(u(abs(u)>pi*B).^2-(pi*B)^2);

end;



SF_recdb=20*log10(abs(SF_rec)/max(abs(SF_rec)));
SF_recdb(SF_recdb<=thres)=thres;

if (taylor_mode==1),
   disp(' ');
   disp(['Scaling factor sigma = ' num2str(sig)]); disp(' ');
   disp('Space Factor nulls located at:');
   disp([num2str(theta_nul','%6.2f') repmat(' deg.',length(theta_nul),1)]);
   HPBW_app=2*asin(sig/(pi*len)*sqrt((acosh(R0)^2-acosh(R0/sqrt(2))^2)))*180/pi;
   
   [y,ind]=min(abs(SF_recdb+3));
   HPBW_graph=2*abs(theta(ind)*180/pi-90);
   disp(' '); disp(['HPBW based on approximate formula = ', num2str(HPBW_app,'%6.2f') ' deg.']);
   disp(['HPBW as computed from the Space Factor graph = ' num2str(HPBW_graph,'%6.2f') ' deg.']);
   
   if output_mode==2,
      if ~exist('fid'),
         fid=fopen(outfile,'wt');
      end;
      fprintf(fid,'%% Taylor synthesis method based on Tschebysheff error\n\n');
      fprintf(fid,['%% Length of line source = ' num2str(len) ' lambda\n']);
      fprintf(fid,['%% Desired constant side lobe level = ' num2str(side_lobe) ' dB\n']);
      fprintf(fid,['%% Desired number of equilevel lobes nbar = ' num2str(nbar) '\n\n']);
      fprintf(fid,['%% Scaling factor sigma = ' num2str(sig) '\n\n']);
      fprintf(fid,'%% Space Factor nulls located at:\n');   
      fprintf(fid,'%% %6.2f deg.\n',theta_nul);
      fprintf(fid,'\n\n%% HPBW based on approximate formula = %6.2f deg.\n',HPBW_app);            
      fprintf(fid,'%% HPBW as computed from the Space Factor graph = %6.2f deg.\n\n',HPBW_graph);
   end;
   
end;

if taylor_mode==2&output_mode==2,
   if ~exist('fid'),
      fid=fopen(outfile,'wt');
      fprintf(fid,'%% Taylor synthesis method based on One-Parameter\n\n');
      fprintf(fid,['%% Length of line source = ' num2str(len) ' lambda\n']);      
      fprintf(fid,['%% Maximum side lobe level = ' num2str(side_lobe) ' dB\n\n']);
   end;   
end;   


% Figure 1
% *********
set(gcf,'units','pixels','position',pos)
plot(zs,abs(Is/max(abs(Is))),'linewidth',2); grid;
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[-len/2, len/2]); 
xlabel('z^{\prime}/\lambda','fontname','timesnewroman','fontsize',18,'verticalalign','top');
ylabel('Current amplitude','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','bottom');

if (taylor_mode==1),
   title(['Current distribution for the Taylor method (Tschebyscheff error, l = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
else
   title(['Current distribution for the Taylor method (One-Parameter, l = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
end;

% Figure 2
% ********
figure(2);
set(gcf,'units','pixels','position',pos);
plot(theta*180/pi,SF_recdb,'linewidth',2); grid;
set(gca,'xlim',[0 180],'fontname','timesnewroman','fontsize',14); 

xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18,'verticalalign','top');
ylabel('Normalized Space Factor (dB)','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','bottom');

h=line([0 180],[-3 -3]); set(h,'color','r','linestyle','--','linewidth',2);
hth=line([0 180],[side_lobe side_lobe],'color',[81 164 109]/255,'linestyle','--','linewidth',2);

set(gca,'units','normalized');
text(0.85,0.95,'-3 dB level','units','normalized','fontname','timesnewroman','fontsize',14);

if (taylor_mode==1),
    title(['Synthesized SF: Taylor method based on Tschebyscheff error (l = ',num2str(len),'\lambda)',], ...
          'fontname','timesnewroman','fontsize',18);
else
   title(['Synthesized SF: Taylor method based on One-Parameter (l = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
end;

% Figure 3
% ********
figure(3);
set(gcf,'units','pixels','position',pos);
elevation([theta theta+pi],[SF_recdb fliplr(SF_recdb)],thres(1),0,4);
set(gca,'fontname','timesnewroman','fontsize',14);
set(findobj(gca,'type','text'),'fontname','timesnewroman','fontsize',14);
set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

if (taylor_mode==1),
   title(['Synthesized SF: Taylor method based on Tschebyscheff error (l = ',num2str(len),'\lambda)',], ...
          'fontname','timesnewroman','fontsize',16)
else   
   title(['Synthesized SF: Taylor method based on One-Parameter (l = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',16);

end;

if ~exist('HPBW_app'),  % one-parameter 
   [y,in]=min(abs(SF_recdb+3));
   HPBW_app=2*abs(90-theta(in)*180/pi);
end;

U=(abs(SF_rec(theta<pi))).^2;
Prad=sum(U.*sin(theta(theta<pi)))*(theta(2)-theta(1));
D=2*max(U)/Prad;
D_mcd=101/(HPBW_app-0.0027*HPBW_app^2);
D_db=10*log10(D);
D_poz=-172.4+191*sqrt(0.818+1/HPBW_app);

if (output_mode==2),
                 
    if taylor_mode==1,
       fidsf=fopen('SF_Taylor.txt','wt');
       fprintf(fidsf,'%% Space factor for the Taylor method based on the Tschebyscheff error\n\n');  
       fprintf(fidsf,'  %% Theta (deg.)      SF (-dB) \n');    
       fprintf(fidsf,'    %6.2f \t      %6.2f\n',[theta*180/pi; SF_recdb]);               
    else
       fidsf=fopen('SF_Taylor.txt','wt');       
       fprintf(fidsf,'%% Space factor for the Taylor method based on One-Parameter\n \n');       
       fprintf(fidsf,'  %% Theta (deg.)      SF (-dB) \n');    
       fprintf(fidsf,'    %6.2f \t      %6.2f\n',[theta*180/pi; SF_recdb]); 
       fprintf(fid,'\n\n%% HPBW based on approximate formula = %6.2f deg.\n\n',HPBW_app);
    end;    
    
    fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Space Factor = %6.2f\n',D);
    fprintf(fid,'%% Directivity (dB) as computed from synthesized Space Factor = %6.2f\n\n',D_db);       
    fprintf(fid,'%% Directivity (dimensionless) from McDonald''s formula = %6.2f\n',D_mcd);
    fprintf(fid,'%% Directivity (dB) from McDonald''s formula = %6.2f\n\n',10*log10(D_mcd));
    fprintf(fid,'%% Directivity (dimensionless) from Pozar''s formula = %6.2f\n',D_poz);      
    fprintf(fid,'%% Directivity (dB) from Pozar''s formula = %6.2f\n',10*log10(D_poz));
        
end;

disp(' ');
if taylor_mode==2,
   disp(['HPBW as computed from the Space Factor graph = ' num2str(HPBW_app,'%6.2f') ' deg.']);
   disp(' ');
end;   
disp(['Directivity (dimensionless) as computed from synthesized Space Factor = ' num2str(D,'%6.2f')]);
disp(['Directivity (dB) as computed from synthesized Space Factor = ' num2str(D_db,'%6.2f')]); 
disp(' ');
disp(['Directivity (dimensionless) from McDonald''s formula = ' num2str(D_mcd,'%6.2f')]);
disp(['Directivity (dB) from McDonald''s formula = ' num2str(10*log10(D_mcd),'%6.2f')]); 
disp(' ');
disp(['Directivity (dimensionless) from Pozar''s formula = ' num2str(D_poz,'%6.2f')]);
disp(['Directivity (dB) from Pozar''s formula = ' num2str(10*log10(D_poz),'%6.2f')]);
disp(' ');  
if output_mode==2,
   disp(['Output data saved in file ''' outfile '''']);
   disp(['Space Factor saved in file ''' pwd '\SF_Taylor.txt''']);
end;   

disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;
fclose('all');



function [y]=fact(x);

for i=1:length(x),
   if (x(i)==0),
      y(i)=1;
   else
      y(i)=prod(1:x(i));
   end;
end;