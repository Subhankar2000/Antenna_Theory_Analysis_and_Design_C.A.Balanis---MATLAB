
%****************************************************************************************************
%   QUARTERWAVE.m
%*************************************************************************
%   This is a MATLAB based program to analyze and design muyltiple-section quarter-wavelength
%   impedance transformers. 
%
%   Quarter-wavelength impedance techniques use either the Binomial or Tschebyscheff method.
%
%   Refer to:  CHAPTER 9 - Quarter-Wavelength Transformers
%
%*******************************************************');
%   THIS PROGRAM IS A MATLAB-BASED SOFTWARE THAT MAKES USE OF THE BINOMIAL AND TSCHEBYSCHEFF DESIGN METHOD
%   OF THE MULTIPLE SECTION QUARTER-WAVELENGTH TRASFORMER FOR IMPEDANCE MATCHING.
%
%   THE GENERAL IDEA OF THE METHOD IS TO INSERT N QUARTER-WAVELENGTH SUBSECTIONS OF TRANSMISSION LINES 
%   BETWEEN THE IMPEDANCE LOAD AND THE CHARACTERISTIC IMPEDANCE OF THE INPUT TRANSMISSION LINE
%   TO PROVIDE AN IMPEDANCE FREQUENCY RESPONSE WITH BINOMIAL/MAXIMALLY-FLAT AND TSCHEBYSCHEFF/EQUAL-RIPPLE 
%   CHARACTERISTICS.
%
%   TO UTILIZE THE METHOD, THE USER MUST DO THE FOLLOWING:
%   1.  SELECT WHICH METHOD TO USE:
%       A.  BINOMIAL
%       B.  TSCHEBYSCHEFF
%   2.  SPECIFY THE:
%       A.  CHARACTERISTIC IMPEDANCE VALUE OF THE INPUT TRANSMISSION LINE
%       B.  LOAD IMPEDANCE; INDICATE WHETHER THE IMPEDANCE IS REAL OR COMPLEX.
%       C.  NUMBER OF QUARTER-WAVELENGTH SUBSECTIONS
%       D.  FRACTIONAL BANDWIDTH 
%               or
%       D'. MAXIMUM TOLERABLE REFLECTION COEFFICIENT (SPECIFY D OR D'; NOT BOTH)  
%
%   THE PROGRAM COMPUTES THE:
%   1.  INTRINSIC REFLECTION COEFFICIENTS AT EACH INTERFACE
%   2.  CHARACTERISTIC IMPEDANCES OF EACH QUARTER-WAVELENGTH SUBSECTION TRANSMISSION LINE 
%   3.  MAXIMUM TOLERABLE REFLECTION COEFFICIENT
%                          or
%   3'. FRACTIONAL BANDWIDTH(COMPUTES 3 OR 3'; NOT BOTH)
%   4.  MAXIMUM INPUT VSWR BASED ON THE MAXIMUM TOLERABLE REFLECTION COEFFICIENT
%
%   THE PROGRAM PLOTS THE MAGNITUDE OF THE INPUT REFLECTION COEFFICIENT AS A FUNCTION 
%   OF THE NORMALIZED FREQUENCY 0<=f/f0<=2.
%
%-----------------------------------------------------------------------------------------------------
%   Written by: PANAYIOTIS I. IOANNIDES - September 2002
%-----------------------------------------------------------------------------------------------------
%
clear all;
close all;
clc; clc;
disp('% ');
disp('%   QUARTER-WAVELENGTH TRANSFORMERS');
disp('%   -------------------------------');
disp('% ');

disp('%   SPECIFY THE METHOD YOU PREFER TO USE');
disp('%   ------------------------------------');
disp('% ');

%selection of the method - Binomial or Tschebyscheff

disp('%   SELECT:');
disp('%   -------');
disp('%   (1) FOR THE BINOMIAL METHOD OR');
disp('%   (2) FOR THE TSCHEBYSCHEFF METHOD'); disp('% ');

method=[]; output_way=[];

while isempty(method)|((method~=1)&(method~=2)),
   method=input('%   THE PREFFERED METHOD IS: ');
end

disp('% ');

disp('%   SELECT OUTPUT MODE:');
disp('%   -------------------');
disp('%   (1) SCREEN');
disp('%   (2) OUTPUT FILE'); disp('% ');  

while isempty(output_way)|((output_way~=1)&(output_way~=2)),
   output_way=input('%   OUTPUT MODE: ');
end

disp('% ');

if (output_way==2),
   filename=[];
   while isempty(filename)
      disp('%   WRITE THE NAME OF THE FILE THAT RESULTS WILL BE WRITTEN IN');
      filename=input('%   (PREFERABLE WITH TXT EXTENSION): ','s');
   end
   disp('% ');
end

%selection of the characteristic impedance of the transmission line

disp('%   GIVE THE CHARACTERISTIC IMPEDANCE OF THE INPUT TRANSMISSION LINE');
disp('%   (IN Ohms, POSITIVE VALUE):');
disp('%   ----------------------------------------------------------------'); disp('% ');
 
 Z0=[];

while isempty(Z0)|(Z0<=0),
   Z0=input('%   CHARACTERISTIC IMPEDANCE OF THE INPUT TRANSMISSION LINE (Ohms): ');
end

%selection of the impedance load

disp('% ');
disp('%   CHOOSE THE TYPE OF THE LOAD IMPEDANCE: ');
disp('%   --------------------------------------'); disp('% ');
disp('%   TYPE (1) REAL');
disp('%   TYPE (2) COMPLEX'); disp('% ');

c=[];

while isempty(c)|((c~=1)&(c~=2)),
   c=input('%   -->: ');
end;

disp('% ');

Rl_rep=[];
Rl=[];
Xl=[];

switch c
case 1
   disp('%   REAL IMPEDANCE LOAD');
   disp('%   -------------------'); disp('% ');
   while isempty(Rl)|(Rl<=0),
     Rl=input('%   GIVE THE VALUE OF THE LOAD RESISTANCE (IN Ohms, POSITIVE VALUE): ');
   end;

otherwise
   disp('%   COMPLEX LOAD IMPEDANCE');
   disp('%   ----------------------'); disp('% ');
    
   while isempty(Rl_rep)|(Rl_rep<0),
      Rl_rep=input('%   GIVE THE VALUE OF THE REAL PART OF THE LOAD (IN Ohms, NOT NEGATIVE VALUE): ');
   end 
  
  disp('% ');
        
   while isempty(Xl)|(Xl==0),
     Xl=input('%   GIVE THE VALUE OF THE IMAGINARY PART OF THE LOAD (IN Ohms, POSITIVE OR NEGATIVE VALUE; NOT 0): ');
   end
   
   result=match(Z0,Rl_rep,Xl);
   s=result(1);
   Rl=result(2);
   
end   

comp=abs((Rl-Z0)/(Rl+Z0));

disp('% '); disp('%   SPECIFY THE NUMBER N OF THE QUARTER-WAVELENGTH SUBSECTIONS (MUST BE INTEGER NUMBER GREATER THAN 1)');
disp('%   --------------------------------------------------------------------------------------------------');  disp('% ');

N=[];

while isempty(N)|(N<=1)|(mod(N,1)~=0),
   N=input('%   (INTEGER NUMBER GREATER THAN ONE) ---> N = ');
end

% selection of either the maximum acceptable bandwidth or the maximum tolerable reflection coefficient

disp('% '); choice=[];
disp('%   SPECIFY EITHER THE FRACTIONAL BANDWIDTH OR');
disp('%   THE MAXIMUM TOLERABLE REFLECTION COEFFICIENT');
disp('%   --------------------------------------------'); disp('% '); 
disp('%   THE FRACTIONAL BANDWIDTH MUST BE');
disp('%   EQUAL OR GREATER THAN 0 AND EQUAL OR LESS THAN 2'); disp('% ');
disp('%   THE MAXIMUM TOLERABLE REFLECTION COEFFICIENT MUST BE');
disp(['%   GREATER THAN 0 AND LESS THAN ',num2str(comp)]); disp('% ');

disp('%   PRESS (1) IF YOU WANT TO SPECIFY THE FRACTIONAL BANDWIDTH'); 
disp('%   PRESS (2) IF YOU WANT TO SPECIFY THE MAXIMUM TOLERABLE REFLECTION COEFFICIENT '); disp('% ');

while isempty (choice)|((choice~=1)&(choice~=2)),
       choice=input('%   ---> choice: ');
end

disp('% ');

switch choice
    case 1
        delta_f=[];
        while isempty(delta_f)|(delta_f<0)|(delta_f>2),
          delta_f= input('%   SPECIFY THE FRACTIONAL BANDWIDTH (0 <= FRACTIONAL BANDWIDTH <= 2): ');
        end
      
     otherwise
        ro_max=[];
        while isempty(ro_max)|(ro_max<=0)|(ro_max>=comp),
          ro_max= input(['%   ENTER THE MAXIMUM TOLERABLE REFLECTION COEFFICIENT (0 < MAXIMUM REFLECTION COEFFICIENT < ',num2str(comp),'): ']);
        end   
 end  
   
   disp('% '); disp('% ');
   
   disp('%   ***********************************************************************');
   disp('%   *-------------------------- R E S U L T S ----------------------------*');
   disp('%   ***********************************************************************');
   
   disp('% '); disp('% ');

%results            

switch method
case 1
   
   disp('%   RESULTS FOR THE BINOMIAL METHOD');
   disp('%   -------------------------------'); disp('% ');
   
otherwise
   
   disp('%   RESULTS FOR THE TSCHEBYSCHEFF METHOD');
   disp('%   ------------------------------------'); disp('% ');

end 

if (output_way==2),
   disp(['%   THE RESULTS ARE ALSO AVAILABLE IN THE FILE ',filename,]);
   disp('%   WHICH IS IN THE SAME FOLDER AS THE MAIN PROGRAM');
end

%output of the characteristic input impedance

disp('% '); disp(['%   CHARACTERISTIC IMPEDANCE OF THE INPUT TRANSMISSION LINE: Z0 = ',num2str(Z0),' Ohms']); disp('% '); 

%output of the load characteristics

switch c
case 1
   disp(['%   THE LOAD IMPEDANCE IS REAL WITH Rl = ',num2str(Rl),' Ohms.']); disp('% ');
otherwise
   if Xl<0,
        disp(['%   THE LOAD IMPEDANCE IS COMPLEX WITH ZL = RL + jXL = ',num2str(Rl_rep),' - j',num2str(abs(Xl)),' Ohms.']); disp('% ');
      else
        disp(['%   THE LOAD IMPEDANCE IS COMPLEX WITH ZL = RL + jXL = ',num2str(Rl_rep),' + j',num2str(Xl),' Ohms.']); disp('% ');
   end
   disp(['%   THE MATCHING WILL BE PERFORMED AT DISTANCE s0 = ',num2str(s)]); 
   disp('%   WAVELENGTHS FROM THE LOAD TOWARDS THE GENERATOR.'); disp('% ');
   disp(['%   THE NEW LOAD IMPEDANCE AT THIS POINT IS REAL: RL = ',num2str(Rl),' Ohms.']); disp('% ');
 end

disp(['%   NUMBER OF SUBSECTIONS USED: ',num2str(N)]); disp('% ');

%results

if (output_way==2),
   fid=fopen(filename,'wt');
   
   switch method
   	case 1
         	fprintf(fid,'\n');
            fprintf(fid,'->   RESULTS FOR THE BINOMIAL METHOD \n');
         	fprintf(fid,'->   ------------------------------- \n \n');
         otherwise
            fprintf(fid,'\n');
         	fprintf(fid,'->   RESULTS FOR THE TSCHEBYSCHEFF METHOD \n');
         	fprintf(fid,'->   ------------------------------------ \n \n');
   	end 
   
   fprintf(fid,'->   CHARACTERISTIC IMPEDANCE OF THE INPUT TRANSMISSION LINE: \n');
   fprintf(fid,['->   Z0 = ',num2str(Z0),' Ohms \n \n']); 
   
   switch c
		case 1
   		fprintf(fid,['->   THE LOAD IMPEDANCE IS REAL WITH Rl = ',num2str(Rl),' Ohms. \n \n']);
   otherwise
         
         if Xl<0,
              fprintf(fid,['->   THE LOAD IMPEDANCE IS COMPLEX WITH ZL = RL + jXL = ',num2str(Rl_rep),' - j',...
              num2str(abs(Xl)),' Ohms. \n \n']);
      	 else
              fprintf(fid,['->   THE LOAD IMPEDANCE IS COMPLEX WITH ZL = RL + jXL = ',num2str(Rl_rep),' + j',...
              num2str(Xl),' Ohms. \n \n']);
         end
         
         fprintf(fid,['->   THE MATCHING WILL BE PERFORMED AT DISTANCE s0 = ',num2str(s),'\n']);
         fprintf(fid,'->   WAVELENGTHS FROM THE LOAD TOWARDS THE GENERATOR. \n \n'); 
         fprintf(fid,['->   THE NEW LOAD IMPEDANCE AT THIS POINT IS REAL: RL = ',num2str(Rl),' Ohms.\n \n']); 
   end

	fprintf(fid,['->   NUMBER OF SUBSECTIONS USED: ',num2str(N),' \n \n']);   

end


switch method
   
   case 1
   
   	switch choice

     		case 1
        	   [dat1 dat2 dat3]=binomial(choice,delta_f,N,Z0,Rl);
        	   disp(['%   THE GIVEN FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1))]); disp('% ');
                   disp(['%   THE ESTIMATED MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: rhom = ',num2str(dat3(2))]); disp('% ');
                  
                   if (output_way==2),
                      fprintf(fid,['->   THE GIVEN FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1)),' \n \n']);
                      fprintf(fid,'->   THE ESTIMATED MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: \n');
                      fprintf(fid,['->   rhom = ',num2str(dat3(2)),' \n \n']);
                   end 
            
                otherwise
                    [dat1 dat2 dat3]=binomial(choice,ro_max,N,Z0,Rl);
                    disp(['%   THE GIVEN MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: rhom = ',num2str(dat3(2))]); disp('% ');
                    disp(['%   THE ESTIMATED FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1))]); disp('% ');
                    if (output_way==2),
                       fprintf(fid,'->   THE GIVEN MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: \n');
                       fprintf(fid,['->    rhom = ',num2str(dat3(2)),' \n \n']); 
                       fprintf(fid,['->   THE ESTIMATED FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1)),' \n \n']); 
                    end             

   	       end

      diff=Rl-dat2(length(dat2));
      
      sign=round(diff/abs(diff));
      
      
      disp('%   THE ESTIMATED INTRINSIC REFLECTION COEFFICIENTS AT EACH INTERFACE');
      disp('%   AND CHARACTERISTIC IMPEDANCES OF EACH QUARTER-WAVELENGTH');
      disp('%   SUBSECTION TRANSMISSION LINE ARE:');
      disp('%   -----------------------------------------------------------------'); disp('% ');
      
      for i=1 : length(dat1),
     		disp(['%   rho',num2str(i-1),' = ',num2str(dat1(i),'%10.6f'),'   Z',num2str(num2str(i-1),'%10.6f'),' = ',...
                num2str(dat2(i),'%10.6f'),' Ohms']);
      end
      
      disp(['%                     Z',num2str(length(dat2)-1),' = ',num2str(dat2((length(dat2))),'%10.6f'),' Ohms']); disp('% ');
      disp(['%   Input VSWR = ',num2str(dat3(3),'%10.6f')]); disp('% ');
      
      if (output_way==2),
         
           fprintf(fid,'->   THE ESTIMATED INTRINSIC REFLECTION COEFFICIENTS AT EACH INTERFACE\n');
	   fprintf(fid,'->   AND CHARACTERISTIC IMPEDANCES OF EACH QUARTER-WAVELENGTH\n');
   	   fprintf(fid,'->   SUBSECTION TRANSMISSION LINE ARE:\n');
   	   fprintf(fid,'->   -----------------------------------------------------------------\n\n');

           for i=1 : length(dat1),
     		   fprintf(fid,['->   rho',num2str(i-1),' =',num2str(dat1(i),'%10.6f'),'   Z'...
                   ,num2str(num2str(i-1),'%10.6f'),' = ',num2str(dat2(i),'%10.6f'),' Ohms \n']);
           end
           fprintf(fid,['->                     Z',num2str(length(dat2)-1),' = ',num2str(dat2((length(dat2))),'%10.6f'),' Ohms \n \n']); 
           fprintf(fid,['->   Input VSWR = ',num2str(dat3(3),'%10.6f \n \n')]); 
         
      end

      switch sign
         case 1
            disp(['%   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms,']);  
            disp(['%   IT IS LESS BY ',num2str(diff),' Ohms BECAUSE THE METHOD IS APPROXIMATE.']);
   	 otherwise
            disp(['%   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms,']);
            disp(['%   IT IS GREATER BY ',num2str(abs(diff)),' Ohms BECAUSE THE METHOD IS APPROXIMATE.']);
      end
      
      if (output_way==2),
         switch sign
         	case 1
               fprintf(fid,['->   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms, \n']);  
               fprintf(fid,['->   IT IS LESS BY ',num2str(diff),' Ohms BECAUSE THE METHOD IS APPROXIMATE. \']);
   	 otherwise
               fprintf(fid,['->   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms, \n']);
               fprintf(fid,['->   IT IS GREATER BY ',num2str(abs(diff)),' Ohms BECAUSE THE METHOD IS APPROXIMATE. \n']);
         end

            fclose(fid);
      end


otherwise,

   
   	switch choice
      	case 1
            [dat1 dat2 dat3]=tschebyscheff(choice,delta_f,N,Z0,Rl);       
            disp(['%   THE GIVEN FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1))]); disp('% ');
            disp(['%   THE ESTIMATED MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: rhom = ',num2str(dat3(2))]); disp('% ');
            if (output_way==2),
                fprintf(fid,['->   THE GIVEN FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1)),' \n \n']);
                fprintf(fid,'->   THE ESTIMATED MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS:\n');
                fprintf(fid,['->   rhom = ',num2str(dat3(2)),' \n \n']);
            end             
      	otherwise
            [dat1 dat2 dat3]=tschebyscheff(choice,ro_max,N,Z0,Rl);  
            disp(['%   THE GIVEN MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: rhom = ',num2str(dat3(2))]); disp('% ');
            disp(['%   THE ESTIMATED FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1))]); disp('% ');
            if (output_way==2),
               fprintf(fid,'->   THE GIVEN MAXIMUM TOLERABLE REFLECTION COEFFICIENT IS: \n');
               fprintf(fid,['->    rhom = ',num2str(dat3(2)),' \n \n']); 
               fprintf(fid,['->   THE ESTIMATED FRACTIONAL BANDWIDTH IS: (Df/f0) = ',num2str(dat3(1)),' \n \n']); 
            end             

   	end
        
      diff=Rl-dat2(length(dat2));
      
      sign=round(diff/abs(diff));
      
      disp('%   THE ESTIMATED INTRINSIC REFLECTION COEFFICIENTS AT EACH INTERFACE');
      disp('%   AND CHARACTERISTIC IMPEDANCES OF EACH QUARTER-WAVELENGTH');
      disp('%   SUBSECTION TRANSMISSION LINE ARE:');
      disp('%   -----------------------------------------------------------------'); disp('% ');

   	for i=1 : length(dat1),
      	    disp(['%   rho',num2str(i-1),' = ',num2str(dat1(i),'%10.6f'),'   Z',num2str(num2str(i-1),'%10.6f'),' = ',...
            num2str(dat2(i),'%10.6f'),' Ohms']);
   	end
   	
      disp(['%                     Z',num2str(length(dat2)-1),' = ',num2str(dat2((length(dat2))),'%10.6f'),' Ohms']); disp('% ');
      disp(['%   Input VSWR = ',num2str(dat3(3),'%10.6f')]); disp('% ');
      
      if (output_way==2),
         
         fprintf(fid,'->   THE ESTIMATED INTRINSIC REFLECTION COEFFICIENTS AT EACH INTERFACE\n');
	      fprintf(fid,'->   AND CHARACTERISTIC IMPEDANCES OF EACH QUARTER-WAVELENGTH \n');
   	   fprintf(fid,'->   SUBSECTION TRANSMISSION LINE ARE: \n');
   	   fprintf(fid,'->   -----------------------------------------------------------------\n \n');

        	for i=1 : length(dat1),
      		fprintf(fid,['->   rho',num2str(i-1),' = ',num2str(dat1(i),'%10.6f'),'   Z',num2str(num2str(i-1),'%8.4f'),' = ',...
                num2str(dat2(i),'%8.4f'),' Ohms \n']);
   		end
         fprintf(fid,['->                     Z',num2str(length(dat2)-1),' = ',num2str(dat2((length(dat2))),'%10.6f'),' Ohms \n \n']);
         fprintf(fid,['->   Input VSWR = ',num2str(dat3(3),'%10.6f \n \n')]);
      end

      
      switch sign
         case 1
            disp(['%   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms,']);
            disp(['%   IT IS LESS BY ',num2str(diff),' Ohms BECAUSE THE METHOD IS APPROXIMATE.']);
         otherwise
            disp(['%   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms,']);
            disp(['%   IT IS GREATER BY ',num2str(abs(diff)),' Ohms BECAUSE THE METHOD IS APPROXIMATE.']);
      end
      
      
      if (output_way==2),
         switch sign
         	case 1
               fprintf(fid,['->   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms, \n']);  
               fprintf(fid,['->   IT IS LESS BY ',num2str(diff),' Ohms BECAUSE THE METHOD IS APPROXIMATE. \n']);
   			otherwise
               fprintf(fid,['->   ALTHOUGH Z',num2str(length(dat2)-1),' SHOULD BE EQUAL TO RL = ',num2str(Rl),' Ohms, \n']);
               fprintf(fid,['->   IT IS GREATER BY ',num2str(abs(diff)),' Ohms BECAUSE THE METHOD IS APPROXIMATE. \n']);
         end       
            fclose(fid);
      end
      
 end 

 disp('% ');
