function [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                    filename,output_mode,fid,poly_opt)


z_nul=z_nul(abs(abs(z_nul)-1)<1e-4);  % use tolerance to get all
						 % nulls on unit circle
psi_nul=angle(z_nul);

disp('Phase of roots on the unit circle (in degrees):')
if ~isempty(psi_nul),
   disp(strjust(num2str(sort(psi_nul*180/pi),'%6.2f'),'right'));
else 
   disp('none');
end;   
disp(' '); 
        	
disp(['Visible region extends from ', num2str(beta*180/pi-360*d) ...
      ' degrees to ', num2str(beta*180/pi+360*d) ' degrees']);
disp(' ');

disp('Roots in invisible region:');
z_invis=z_nul(psi_nul<(beta-2*pi*d)|psi_nul>(beta+2*pi*d));

if isempty(z_invis),
   disp('none'); disp(' ');
else
   disp(char(cosm(z_invis))); disp(' ');
end; 

disp('Roots in visible region:');            
z_visi=z_nul(psi_nul>(beta-2*pi*d)&psi_nul<(beta+2*pi*d));
           
if isempty(z_visi),              
   disp('none'); disp(' ');
else            
  disp(char(cosm(z_visi))); disp(' ');           
end;

if (output_mode==2),             
    fprintf(fid,'\n');              
    fprintf(fid,'%% Phase of roots on the unit circle (in degrees):\n');
    if ~isempty(psi_nul),
       fprintf(fid,'%% %-+6.2f\n',psi_nul*180/pi);              
    else
       fprintf(fid,'%% none\n');
    end;   
    fprintf(fid,'\n');              
    fprintf(fid,['%% Visible region extends from %6.2f degrees to ' ...
                 '%6.2f degrees\n'], beta*180/pi-360*d,beta*180/pi+360*d);
    fprintf(fid,'\n');
    fprintf(fid,'%% Roots in invisible region:\n'); 
              
    if isempty(z_invis),         
        fprintf(fid,'%% none\n');                  
    else
        fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_invis.'); imag(z_invis.')]);
    end;
                           
    fprintf(fid,'\n%% Roots in visible region:\n');

    if isempty(z_visi),	         
        fprintf(fid,'%% none\n');
    else
        fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_visi.'); imag(z_visi.')]);   
    end;

end;

psi_nul_vis=psi_nul((abs(psi_nul-beta)<=(2*pi*d)+1e-3));   
                       % psi's in visible region (again use tolerance)
	   
if (isempty(psi_nul_vis)),
    disp('No nulls exist');

    if (output_mode==2),
        fprintf(fid,'\n%% No nulls exist \n');
    end;
	   
else

% get only the real part of nulls 	       
theta_nul=real(acos((psi_nul_vis-beta)/(2*pi*d)));

    for mu=1:ceil(d*(1-min(cos(theta_nul)))), % extra code to cover 
                                              % the case d>0.5 lambda  
                 
        theta_nul1(mu,:)=acos(cos(theta_nul')+mu/d);
        theta_nul2(mu,:)=acos(cos(theta_nul')-mu/d);
    end;

    theta_nul=[theta_nul; theta_nul1(:); theta_nul2(:)];
    theta_nul(abs(imag(theta_nul))>=1e-4)=[];
    nul_tot=[theta_nul*180/pi];

	       
    for w=1:length(nul_tot),   % extra code to remove multiple null
		           % occurences
       nul_final(w)=nul_tot(1);
       nul_tot(abs(nul_tot-nul_final(w))<=1e-2)=[];
	          
       if (isempty(nul_tot)),
          break;
       end;
    end;	

	       
    if (isempty(nul_final)),
       disp('No Array Factor nulls exist');
    else
       ats=nul_final(find(nul_final(nul_final<0|nul_final>180)));
       disp(ats);
    end;  
       
    if ~isempty(ats),
        disp('Null found in invisible region');
    end;

    disp('Array Factor nulls inside visible region were found at:');
    disp(strcat(num2str(sort(nul_final'),'%6.2f'),  ...
         repmat(' deg.',length(nul_final),1)));
    
    if exist('fid')&(output_mode==2),  
       fprintf(fid,'\n%% Array Factor nulls inside visible region were found at:\n');  
       fprintf(fid,'%% %6.2f deg.\n',sort(nul_final'));
    end;
             
end;






