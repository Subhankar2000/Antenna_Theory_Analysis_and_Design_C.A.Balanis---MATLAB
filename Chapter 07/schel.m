% this program synthesizes a linear array using the Schelkunoff
% polynomial method

% Written by Marios Gkatzianas	
% Arizona State University, September 2002

function []=schel(output_mode,filename,pos);

% load bal.mat;
[x,map]=imread('bal.tiff');

close all;
warning off;

mode_schel=3;
while (mode_schel~=1)&(mode_schel~=2),
    mode_schel=input(['Schelkunoff synthesis options\n', ...
           '*****************************\n',... 
           '1. Compute polynomial for specific nulls\n', ...
           '2. Compute nulls for a specified polynomial\n', ...
           'Choose one of the above:\n','->']);
     
    if isempty(mode_schel)|((mode_schel~=1)&(mode_schel~=2)),
        hso=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal'); 
    end;   
end;

d=0;
while ~isreal(d)|(isempty(d)|(d<=0)),
   d=input('Specify spacing between elements (in wavelengths)\n');
   d=d+1e-12;
   if ~isreal(d)|isempty(d)|d<=0,
      hss=msgbox('Spacing must be positive','Invalid choice','custom',x,map,'modal');
   end;   
end;

if (mode_schel==1),
   disp('It is assumed that beta=0'); beta=0;
else

   beta=[];
   while ~isreal(beta)|isempty(beta),	
      beta=input('Specify progressive phase shift beta (in degrees)\n');
      if ~isreal(beta)|isempty(beta),
         create.WindowStyle='modal';
         create.Interpreter='tex';
         hsp=msgbox(strvcat('Progressive phase \beta must be non-empty','and real'), ...
                           'Invalid choice','custom',x,map,create);
      end;   
   end;

end;
beta=beta*pi/180;
poly_opt=4;
switch mode_schel
    case 1,  % nulls -> poly

        Nnul=0;
        while ~isreal(Nnul)|isempty(Nnul)|Nnul<=0|mod(real(Nnul),1)~=0,
           Nnul=input('Specify number of desired nulls\n');
           if ~isreal(Nnul)|isempty(Nnul)|isempty(Nnul)|Nnul<=0|mod(real(Nnul),1)~=0,
              hnnul=msgbox('Number of nulls must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');
           end;   
        end;

	Nel=Nnul+1;
	nul=zeros(size(Nnul));

	for n=1:Nnul,
	    aux=[];

	    while ~isreal(aux)|isempty(aux)|aux<0|aux>180,
     	     aux=input(['Specify position of null #',num2str(n), ...
                   ' (in degrees)\n']);
           if ~isreal(aux)|isempty(aux)|aux<0|aux>180,
              hposn=msgbox('Position of null must be between 0 and 180 degrees','Invalid choice', ...
                            'custom',x,map,'modal'); 
           end;     
       end;

	    nul(n)=aux;
	end;

	nul=nul*pi/180;   % convert to rads

    case 2, % poly -> nulls   
     

	while ((poly_opt~=1)&(poly_opt~=2)&(poly_opt~=3)),
	    poly_opt=input([ 'Polynomial definition options\n', ...         
                 '*****************************\n', ...
                 '1. Enter vector of polynomial coefficients in descending form\n', ...
                 '---------------------------------------------------\n', ...
                 'Example: z^3-z^2+z-1 should be input as [1 -1 1 -1]\n', ...
                 '---------------------------------------------------\n\n', ...
                 '2. Specify polynomial as a function of z\n', ... 
                 '--------------------\n', ...
                 'Example: z^3-z^2+z-1\n', ...
                 '--------------------\n\n', ...
                 '3. Specify roots of the desired polynomial\n' ...
                 '--------------------\n', ...
                 'Example(4 roots): +j, +1, -j, 0.707 + 0.707j\n',... 
                 '--------------------\n', ...
                 'Choose one of the above\n','->' ]);
           
       if isempty(poly_opt)|((poly_opt~=1)&(poly_opt~=2)&(poly_opt~=3)),
          hpoly=msgbox('Specified number must be either 1,2 or 3','Invalid choice', ...
                        'custom',x,map,'modal');
       end;   
   end;

	if (poly_opt==1),
	   a=[];
	   while isempty(a),
	      a=input(strcat('Enter row vector of polynomial coefficients', ...
                        ' (in descending powers of z)\n','Example: z^3-z^2+z-1 should be input as [1 -1 1 -1]\n'));	
         if isempty(a)|~isa(a,'double'),
            hrow=msgbox(strvcat('Polynomial coefficients must be inside square brackets'), ...
                        'Invalid choice', 'custom',x,map,'modal');
         end;   
	   end;

	   a=fliplr(a);   % flip a in ascending powers of z   
       
   elseif (poly_opt==2),
       pstring=[];
       while isempty(pstring)|~isempty(find(pstring(isletter(pstring))~='z'&pstring(isletter(pstring))~='i')),
          
           pstring=input(strcat('Specify polynomial as a function of z. ', ...
              ' Use ^ for all powers and i for the imaginary unit and don''t put quotes\n', ...
              'Example: z^3-z^2+z-1\n'),'s');
           if isempty(pstring)|~isempty(find(pstring(isletter(pstring))~='z'&pstring(isletter(pstring))~='i')),                                                      
              hpstr=msgbox(strvcat('String must be non empty and must only contain','a z variable (apart from i)'), ...
                           'Invalid choice','custom',x,map,'modal');             
           end;   
       end;
       pz=inline(vectorize(pstring));

   else  % define roots

       Nroot=[];   
       while (isempty(Nroot)|Nroot<=0|~isreal(Nroot)|mod(Nroot,1)~=0),    
          Nroot=input('Specify number of roots for the polynomial\n');
          if isempty(Nroot)|((Nroot<=0)|~isreal(Nroot)|mod(Nroot,1)~=0),
              hnumroot=msgbox('Number or roots must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');   
          end;
       end;

       root=[];

       for n=1:Nroot,
           temp_root=[];               
           while isempty(temp_root),
              temp_root=input(['Specify root #', num2str(n),' of the polynomial\n']);
              
              if isempty(temp_root),
                 hroot=msgbox(['Root #' num2str(n) ' must be non-empty'],'Invalid choice', ...
                               'custom',x,map,'modal');
              end;   
           end;

           root=[root temp_root];
           
% the following portion of the code includes the conjugate roots to make the polynomial real. 
% To be used optionally           
           
%           if imag(temp_root)~=0,
%               disp(['The root ', num2str(conj(temp_root)) ' is also included to' ...
%                     ' make the polynomial have real coefficients']);
%               root=[root conj(temp_root)];
%           end;
        end;
                
        polycoef=poly(root);
        a=polycoef;
        
        spoly=[];
        for n=1:length(a)-1,
           stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' '*z^' num2str(length(a)-n) ' '];
           spoly=[spoly stemp];
        end;
        
        spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];
                
    end;

end;  % end of switch 

Ntheta=0;

while ~isreal(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
   Ntheta=input(['Specify number of sample points for the angle theta', ...
         ' for the pattern plots\n', '[ENTER for default=360]: ']); 
   Ntheta(isempty(Ntheta))=360;
   
   if ~isreal(Ntheta)|Ntheta<=0|mod(Ntheta,1)~=0,
      create.WindowStyle='modal';
      create.Interpreter='tex';
      hNthe=msgbox('Number of sample points for \theta must be integer positive','Invalid choice', ...
         'custom',x,map,create);      
   end;   
end;

thres=0;
while ~isreal(thres)|thres>=0,
   thres=input(['Specify lowest threshold of Array Factor for plotting' ...
         ' purposes (in -dB)\n','[ENTER for default=-120]:']);
   thres(isempty(thres))=-120;
   
   if ~isreal(thres)|thres>=0,
      hthres=msgbox('Threshold must be negative','Invalid choice','custom',x,map,'modal');
   end;
   
end;
thres(isempty(thres))=-120;

theta=linspace(0,2*pi,2*Ntheta);
psi=2*pi*d*cos(theta)+beta;


switch mode_schel,
    case 1   % nulls -> polynomial	

	     psi_nul=2*pi*d*cos(nul)+beta;
	     a=poly(exp(j*psi_nul));  % coefficients of z polynomial
        
        spoly=[];
        for n=1:length(a)-1,
           stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' '*z^' num2str(length(a)-n) ' '];
           spoly=[spoly stemp];
        end;
        
        spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];       
        
	     z=exp(j*psi);

	     for n=1:Nel,
	        temp(n,:)=z.^(n-1);
	     end;

	     AF=fliplr(a)*temp;  % final array factor
        	
        for mu=1:ceil(d*(1-min(cos(nul)))),    % extra code to cover the case
					       % d>0.5 lambda  
           theta_nul1(mu,:)=acos(cos(nul)+mu/d)*180/pi;
           theta_nul2(mu,:)=acos(cos(nul)-mu/d)*180/pi;
        end;

	     theta_nul=[theta_nul1(:) theta_nul2(:)];
        theta_nul(abs(imag(theta_nul))>=1e-4)=[];
	     nul_tot=[nul*180/pi theta_nul];

	     for w=1:length(nul_tot),   % extra code to remove multiple null
				   % occurences
	         nul_final(w)=nul_tot(1);
            nul_tot(abs(nul_tot-nul_final(w))<=1e-2)=[];
	         if (isempty(nul_tot)),
                break;
	         end;
	     end;
        
        disp(' '); disp('Schelkunoff polynomial found as:');
        disp(spoly); disp(' ');
	     disp('Array Factor nulls found at:');
		  disp(strcat(num2str(sort(nul_final)','%6.2f'), ...
             repmat(' deg.',length(nul_final),1)))                  
        fprintf(2,'\nExcitation coefficients:\n');   % save output on screen too
        fprintf(2,'   Real part \t    Imag part\n');
        fprintf(2,'%10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
  
          
        if output_mode==2,
           if ~exist('fid'),
              fid=fopen(filename,'wt');
           end;
           fprintf(fid,'%% Schelkunoff synthesis method\n\n');
           fprintf(fid,['%% Spacing d=' num2str(d) ' lambda\n']);
           fprintf(fid,['%% Progressive phase beta=' num2str(beta*180/pi) ' degrees\n']);             
           fprintf(fid,'\n%% User specified nulls:\n');
           fprintf(fid,'%% %5.2f\n',nul'*180/pi);
           fprintf(fid,'\n%% Schelkunoff polynomial:\n');
           fprintf(fid,'%% %s\n',spoly);  
           fprintf(fid,'\n%% Array Factor nulls found at:\n');
           fprintf(fid,'%% %5.2f\n',sort(nul_final'));
           fprintf(fid,'\n%% Excitation coefficients:\n');
           fprintf(fid,'%% Real part \t Imag part\n');
           fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);                                  
        end;

    case 2  % polynomial -> nulls

        if (poly_opt==3),
            psi=2*pi*d*cos(theta)+beta;
            z=exp(j*psi);
            AF=polyval(polycoef,z);
            strpoly=[];
            z_nul=root.';
            Nel=length(z_nul)+1;            

            for n=1:length(polycoef)-1,
               if polycoef(n)~=0
                  tempstr=[num2str(polycoef(n),'%+5.2f') '*z^' ...
                           num2str(length(polycoef)-n)];
               end; 
  	            strpoly=[strpoly tempstr];                
            end;
            
            if polycoef(end)~=0, 
               strpoly=[strpoly num2str(polycoef(end),'%+5.2f')];
            end;  

            disp('Roots of the polynomial:');
            disp(char(cosm(root.')));
            disp('The Schelkunoff polynomial is:');           
            disp(spoly);

	         if (output_mode==2),
                fid=fopen(filename,'wt');
                fprintf(fid,'%% Schelkunoff synthesis method\n');
                fprintf(fid,['%% Spacing d=' num2str(d) ' lambda\n']);
                fprintf(fid,['%% Progressive phase beta=' num2str(beta*180/pi) ' degrees\n\n']);
                fprintf(fid,'%% Roots of the polynomial:\n');
                fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_nul.'); imag(z_nul.')]);         
                fprintf(fid,'\n%% The Schelkunoff polynomial is:\n');
                fprintf(fid,'%% %s\n',spoly);
                fprintf(fid,'\n%% Excitation coefficients:\n');
                fprintf(fid,'%% Real part \t Imag part\n');
                fprintf(fid,'%10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
            end;      

            if ~exist('fid'),
	             fid=[];          
            end;

            [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                     filename,output_mode,fid,poly_opt);

        end;  % poly_opt=3

	     if (poly_opt==1),  % row vector
	        z=exp(j*psi);
           Nel=length(a);

	        for n=1:Nel,
	           temp(n,:)=z.^(n-1);
           end;                      
           
           AF=fliplr(a)*temp;
 
	        z_nul=roots(a);
           a=fliplr(a);
           
           spoly=[];
           for n=1:length(a)-1,
              stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' ...
                     '*z^' num2str(length(a)-n) ' '];
               spoly=[spoly stemp];
           end;
        
           spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];
           
           disp(' ');
           disp('Schelkunoff polynomial is:');
           disp(spoly); disp(' ');                                                  
           
           disp('Roots of the polynomial:');
           disp(char(cosm(z_nul)));                      
                                          
           if (output_mode==2),
              fid=fopen(filename,'wt');
              fprintf(fid,'%% Schelkunoff synthesis method\n\n');
              fprintf(fid,['%% Spacing d= ' num2str(d) ' lambda\n']);
              fprintf(fid,['%% Progressive phase beta= ' num2str(beta*180/pi) ' degrees\n\n']);                      
           
              strpoly=[];
              for n=length(a):-1:2,
                 if a(n)~=0
                    tempstr=[num2str(a(n),'%+5.2f') '*z^' ...
                             num2str(n-1)];
                 end; 
  	              strpoly=[strpoly tempstr];                
              end;
           
              if a(1)~=0, 
                 strpoly=[strpoly num2str(a(1),'+%6.2f')];
              end; 
                     
              fprintf(fid,'%% The Schelkunoff polynomial is:\n%% %s\n\n',spoly);
              fprintf(fid,'%% Roots of the polynomial:\n');
              fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_nul.'); imag(z_nul.')]);  
              fprintf(fid,'\n%% Excitation coefficients:\n');
              fprintf(fid,'%% Real part \t Imag part\n');
			     fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);        
           end;                     

           if ~exist('fid'),
	          fid=[];          
           end;

           [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d,filename,output_mode,fid,poly_opt);                                 
                 
        elseif (poly_opt==2),   % inline function
                                % attention:  add more code for output file       
                                                                        
	        z=exp(j*psi);
           AF=pz(z);
           
           strpz=char(pz);
           strpz(strpz=='.')=[]; 
           
           if exist('sym')~=2,
              disp(['Symbolic toolbox not installed. Skipping symbolic null computations.' ...
                    ' You have to determine the nulls visually.']);
              disp(' '); disp('Schelkunoff polynomial is:');
              disp(strpz); disp(' ');
              
           else   
              zs=sym('zs');                                          
              pzs=strrep(char(pz),'z','zs');
              pzs(pzs=='.')=[];
              zpoly=sym2poly(sym(pzs));
              z_nul=roots(zpoly);           
              a=zpoly;                                         
                                         
              disp(' '); disp('Schelkunoff polynomial is:');
              disp(strpz); disp(' ');            
              disp('Excitation coefficients:');
              disp(char(cosm(num2str(a.','%8.4f')))); disp(' ');
                            
              if ~exist('fid')&(output_mode==2),
                 fid=fopen(filename,'wt');
                 fprintf(fid,'%% Schelkunoff synthesis method\n\n');
                 fprintf(fid,['%% Spacing d= ' num2str(d) ' lambda\n']);
                 fprintf(fid,['%% Progressive phase beta= ' num2str(beta*180/pi) ' degrees\n\n']);                 
                 fprintf(fid,'%% The Schelkunoff polynomial is:\n%% %s\n',strpz);
                 fprintf(fid,'\n%% Excitation coefficients:\n');
                 fprintf(fid,'%% Real part \t Imag part\n');
                 fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
              end;
              
              if ~exist('fid'),
                 fid=[];
              end;
              
              [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                     filename,output_mode,fid,poly_opt);                                 
              Nel=length(a);
           end;                     
           
	     end;   % if 

end;

	AF_db=20*log10(abs(AF)/max(abs(AF)));
	AF_db(AF_db<=thres)=thres;
   
% directivity computation
% ***********************   

dth=theta(2)-theta(1);   % theta in radians
U=(abs(AF(theta<pi))).^2;
Prad=2*pi*sum(U.*sin(theta(theta<pi))*dth);
D=4*pi*max(U)/Prad;
D_db=10*log10(D);
   
%if (mode_schel==1),
if exist('a'),

      % Figure 1
      % ********
      figure(1);
      set(gcf,'units','pixels','position',pos);
      subplot(2,1,1); stem(fliplr(abs(a))); grid on; 
      set(gca,'fontname','timesnewroman','fontsize',14,'xtick',[1:Nel]);	
      ylabel('Amplitude excitation','fontname','timesnewroman','fontsize',18);
      title(['Amplitude and phase excitation (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontsize',18);        
      set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

      subplot(2,1,2); stem(fliplr(angle(a))*180/pi); grid on;
      set(gca,'fontname','timesnewroman','fontsize',14);   
      set(gca,'xtick',[1:Nel]);	set(gca,'ylim',[-180 180]);
%      set(gca,'ytick',[-180 -90 0 90 180]);
      xlabel('element number','fontname','timesnewroman','fontsize',18);
      ylabel('Phase excitation (in degrees)','fontname','timesnewroman','fontsize',18);  
      set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

end;

% Figure 2
% ********
figure(2);
set(gcf,'units','pixels','position',pos);
set(gca,'fontname','timesnewroman','fontsize',14); 
hl=plot(theta*180/pi,AF_db); grid; zoom; thres=get(gca,'ylim'); grid on;
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
set(hl,'linewidth',2); 
xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18, ...
       'verticalalignment','top');
ylabel('Normalized Array Factor (in dB)','fontname','timesnewroman', ...
       'fontsize',18,'verticalalignment','bottom');
title(['Synthesized Array Factor using Schelkunoff polynomial(N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
      'fontname','timesnewroman','fontsize',18);


% Figure 3
% ********
figure(3);

if (mod(thres(1),4)~=0),  
   thres(1)=ceil(thres(1)/4)*4;   
end;

set(gcf,'units','pixels','position',pos);

elevation(theta,AF_db,thres(1),0,4); 
set(gca,'fontname','timesnewroman','fontsize',14);
set(findobj(gca,'type','text'),'fontname','timesnewroman','fontsize',14);
set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

title(['Synthesized Array Factor using Schelkunoff polynomial(N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
      'fontname','timesnewroman','fontsize',18);
   
disp(' ');
disp(['Directivity (dimensionless) as computed from synthesized Array Factor = ' num2str(D,'%6.2f')]);
disp(['Directivity (dB) as computed from synthesized Array Factor = ' num2str(D_db,'%6.2f\n')]);
           
if (output_mode==2),

   if exist('poly_opt')&poly_opt==1&~isempty(nul_final),
       fprintf(fid,'\n%% Array Factor nulls found at:\n');
       fprintf(fid,'%% %6.2f deg. \n',sort(nul_final'));
   end;

   if ~exist('fidm'),
      fidm=fopen('AF_array.txt','wt');
   end;
   
%% WARNING: add more code to write relevant information in second file   
  

   fprintf(fidm,'%% Schelkunoff polynomial:\n');
   if poly_opt==2,
      fprintf(fidm,['%% ' strpz '\n']);
   elseif mode_schel==1|poly_opt==1|poly_opt==3,
      fprintf(fidm,['%% ' spoly '\n']);
   end;      
   fprintf(fidm,'\n%s\n','% Theta (deg.)    AF (dB)');    
   fprintf(fidm,'   %6.2f \t %6.2f\n', [(theta*180/pi); AF_db]);
   fprintf(fidm,'\n');   
   fclose(fidm);
   fprintf(fid,'\n%% Directivity (dimensionless) as computed from Array Factor = %6.2f\n',D);
   fprintf(fid,'%% Directivity (dB) as computed from Array Factor = %6.2f\n',D_db);  
   fclose(fid);
   
   disp(['Output saved in file ''',filename,'''.']);
   disp(['Array Factor saved in file ''' pwd '\AF_array.txt''']);   
end;  

disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;

warning on;