function Gamma_comp=tschebcomp(choice,p,N,Z0,Rl);

% QUARTER WAVELENGTH TRANSFORMER WITH THE USE OF THE TSCHEBYSCHEFF METHOD
%
% INPUT PARAMETERS ARE EITHER THE MAXIMUM FRACTIONAL BANDWIDTH OR THE 
%
% MAXIMUM TOLERABLE REFLECTION COEFFICIENT, THE NUMBER OF USED
%
% SUBSECTIONS (N), THE CHARACTERISRIC INPUT IMPEDANCE (Z0) AND
%
% THE LOAD INPUT RESISTANCE (Rl)

mult1=(Rl-Z0)/(Rl+Z0);

%derivation of the Tschebyscheff polynomials up to the Nth order

T=zeros(N+1);
	T(1,N+1)=1;
	T(2,N)=1;

	ord(1:N)=2:N+1;

	ord(N+1)=1;

	for i = 2 : N,
   	T(i+1,:)=2*T(i,ord)-T(i-1,:);
	end
   
TN=T(N+1,:); %TN is the Nth order Tschebyscheff polynomial
   
switch choice
case 1
   
   	%estimation of the maximum reflection coefficient if the fractional bandwidth is known
      
      thetam=(pi/2)*(1-p/2);
   	secthetam=1/cos(thetam);
   	dat3(1)=p;
      dat3(2)=abs(mult1/polyval(TN,secthetam));
      
   otherwise
      
      %estimation of the fractional bandwidth if the maximum reflection coefficient is known
      
		polynom=T(N+1,:);
	
		zero_power=abs((1/p)*mult1);

		sub=T(N+1,N+1)-zero_power;
	
		polynom(N+1)=sub;

		rts=roots(polynom); %solution of the Nth order polynomial for a certain value
	
		pos=find((real(rts)>0)&(imag(rts)==0));
	
      secthetam=min(rts(pos));
   
      secthetam2=cos((1/N)*acos(abs((1/p)*mult1)));   
      
      costhetam=1/secthetam;
   
      thetam=acos(costhetam);
   
      dat3(1)=2*(1-(2/pi)*thetam);
      
      dat3(2)=p;

end

f=0 : 0.001 : 2; %range of the fractional bandwidth

theta=(pi/2)*f;

mult2=(mult1)/(polyval(TN,secthetam));

exponent=exp(-j*N*theta);

mult3=polyval(TN,(secthetam*cos(theta)));

Gamma_comp=abs(mult2.*exponent.*mult3); %estimation of the Input reflection coeffecient in all the arange of fractional bandwidth

for i = 1: 1 : N+1,
   w(i)=secthetam^(N+1-i);
end  

wTN=(mult2/2)*w.*TN;

rt_test=mod(N,2); %determination if N is even or odd

ii=0; jj=0;

switch rt_test,
   
case 0 %N is even
   rt=N/2+1;
     
   for ii=1 : rt,  
      for jj=1:rt,
         MAT(ii,jj)=T(N+3-2*ii,2*jj-1);
      end
   end

   MAT=MAT';  
   
   ind=find(wTN~=0);
   
   new_wTN=(wTN(ind))';
   
   roots_tsch=(MAT^(-1)*new_wTN)';
   
   roots_tsch(rt)=2*roots_tsch(rt);

   ro=[roots_tsch fliplr(roots_tsch(1:((length(roots_tsch))-1)))]; %estimation of the reflections coefficients

otherwise %N is odd
   
   rt=(N+1)/2;
     
   for ii=1 : rt,  
      for jj=1:rt,
         MAT(ii,jj)=T(N+3-2*ii,2*jj-1);
      end
   end

   MAT=MAT';  
   
   ind=find(wTN~=0);
   
   new_wTN=(wTN(ind))';
   
   roots_tsch=(MAT^(-1)*new_wTN)'; 
   
   ro=[roots_tsch fliplr(roots_tsch)]; %estimation of the reflections coefficients

   
 end
 
 dat1=ro;
 
 Z=zeros(1,N+2);

Z(1)=Z0;

Z(N+2)=Rl;

for i=2 : N+2,
   
   Z(i)=((1+ro(i-1))./(1-ro(i-1)))*Z(i-1); %estimation of the characteristic impedance of each subsection
   
end;

dat2=Z;

dat3(3)=(1+ dat3(2))/(1- dat3(2)); %estimation of the input VSWR

vert_plot=0 : 0.01 : abs(mult1);

max_plot(1:length(f))=dat3(2); clear j;

Gamma_single=abs(.5*mult1*(1+exp(-j*2*theta)));
