function y = U(THETA, PHI)
%  **Example**
%y = sin(THETA).*(sin(PHI)).^2;
%y = (cos(THETA)).^2.*(sin(PHI)).^2+(cos(PHI)).^2;
%y = (1-(sin(THETA)).^2.*(sin(PHI)).^2);
%y = sin(THETA).^3;
%   ****Example****
%PI=3.1415927;
%k=2*PI;
%d=1/2;
%BETA=-k*d;
%PSI=k*d*cos(THETA+BETA);
%N=4;
%y =(sin(N*PSI/2)./(N*sin(PSI/2))).^2; 
%y =(sin(N*PSI/2)./(N*(PSI/2))).^2;
%   *****Example********
%a1=4.9886;
%a2=3.0851;
%a3=1;
%d=0.7124;
%u=pi*d*cos(THETA);
%AF=a1*cos(u)+a2*cos(3*u)+a3*cos(5*u);
%E=cos(pi*cos(THETA)/2)./sin(THETA);
%y=abs(E.*AF).^2;
%y=abs(AF).^2;
N=5;
p=1.0;
LO=1.03061;
S=0.2493;
PI=3.1415927;
PSI=2*PI*(S*cos(THETA)-LO/p);
y=(cos(THETA).*(sin(N*PSI/2)./sin(PSI/2))).^2;



%