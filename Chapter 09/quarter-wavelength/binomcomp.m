function Gamma_comp=binomcomp(choice,p,N,Z0,Rl);

% QUARTER WAVELENGTH TRANSFORMER WITH THE USE OF BINOMIAL METHOD
%
% INPUT PARAMETERS ARE EITHER THE MAXIMUM FRACTIONAL BANDWIDTH OR THE 
%
% MAXIMUM TOLERABLE REFLECTION COEFFICIENT, THE NUMBER OF USED
%
% SUBSECTIONS (N), THE CHARACTERISRIC INPUT IMPEDANCE (Z0) AND
%
% THE LOAD INPUT RESISTANCE (Rl)

mult1=(Rl-Z0)/(Rl+Z0);

n=0 : N;

C_n=(fact(N))./(fact(N-n).*fact(n));

ro_n=2.^(-N)*((Rl-Z0)/(Rl+Z0)).*C_n; %estimation of the reflections coefficients

Z=zeros(1,N+2);

Z(1)=Z0;

Z(N+2)=Rl;

for i=2 : N+2,
   
   Z(i)=((1+ro_n(i-1))./(1-ro_n(i-1)))*Z(i-1); %estimation of the characteristic impedance of each subsection
   
end;

dat1=ro_n; 

dat2=Z;

switch choice
case 1
   dat3(1)=p;
   dat3(2)=abs(((cos((pi/4)*(2-p)))^N)*(mult1));
otherwise
   dat3(2)=p;
   dat3(1)=abs(2-(4/pi)*acos((p/abs(mult1))^(1/N)));
end

dat3(3)=(1+ dat3(2))/(1- dat3(2)); %estimation of the input VSWR

f=0 : 0.001 : 2;

theta=(pi/2)*(f);

exptheta=repmat(theta,[length(n) 1]);

expn=repmat(n',[1 length(f)]);

coeff=(exp(-j.*2*expn.*exptheta));

expro_n=repmat((ro_n)',[1 length(f)]);

Gamma_comp=abs(sum(expro_n.*coeff,1)); %estimation of the Input reflection coeffecient in all the arange of fractional bandwidth

vert_plot=0 : 0.01 : abs(mult1);

max_plot(1:length(f))=dat3(2); clear j;

Gamma_single=abs(.5*mult1*(1+exp(-j*2*theta)));
