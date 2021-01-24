function [theta,phi]=est_doa(X,M,N,dx,dy,num_sig)

%*************************************************************************
%  This function computes the directions of arrival of all signals
%  given the array data, the number of elements, element spacing and
%  number of signals. The function works for a uniform planar or uniform linear array.  
%  If the array is planar, the 2-D Unitary Estimation of Signal Parameters via
%  Rotational Invariance Techniques(ESPRIT) algorithm is used.  If the array 
%  is linear, the total least squares(TLS) ESPRIT is used.

%  Input Variables
%  X        -  MNxL data matrix where L is the number of samples per element.
%              Each row of X contains samples from 1 element.
%  M        -  Number of elememts in x direction.
%  N        -  Number of elements in y direction(for linear array either M=1 or N=1.)
%  dx       -  Spacing between elements in wavelengths in x direction.
%              For linear array, dx is used as interelement spacing and dy is ignored.
%  dy       -  Spacing between elements in wavelengths in y direction.
%  num_sig  -  Number of signals present.

%  Output variables
%  theta    -  Vector containing azimuth direcions of arrival for all signals in radians.
%  phi      -  Vector containing elevation directions of arrival for all signals in radians.
%              If array is linear then phi will be an empty variable.
%*************************************************************************


if N~=1 & M~=1  %if N=1 or M=1 then the array is linear, if not, then use Unitary ESPRIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations for Unitary ESPRIT for planar array
% See Zoltowski, Haardt and Mathews'  paper(Feb. 1996 IEEE 
% Trans. on Sig. Proc.) for a derivation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up constants used in Unitary ESPRIT
if rem(N,2)==0  %check if N is even
QN=[eye(N/2) j*eye(N/2) ; fliplr(eye(N/2)) -j*fliplr(eye(N/2))];
QNM1=[eye(N/2-1) zeros(N/2-1,1) j*eye(N/2-1) ; zeros(N/2-1,1)' sqrt(2) zeros(N/2-1,1)' ; fliplr(eye(N/2-1)) zeros(N/2-1,1) -j*fliplr(eye(N/2-1)) ];
else
  QN=[eye((N-1)/2) zeros((N-1)/2,1) j*eye((N-1)/2) ; zeros((N-1)/2,1)' sqrt(2) zeros((N-1)/2,1)' ; fliplr(eye((N-1)/2)) zeros((N-1)/2,1) -j*fliplr(eye((N-1)/2)) ];
  QNM1=[eye((N-1)/2) j*eye((N-1)/2) ; fliplr(eye((N-1)/2)) -j*fliplr(eye((N-1)/2))];
end

if rem(M,2)==0  %check if M is even
QM=[eye(M/2) j*eye(M/2) ; fliplr(eye(M/2)) -j*fliplr(eye(M/2))];
QMM1=[eye(M/2-1) zeros(M/2-1,1) j*eye(M/2-1) ; zeros(M/2-1,1)' sqrt(2) zeros(M/2-1,1)' ; fliplr(eye(M/2-1)) zeros(M/2-1,1) -j*fliplr(eye(M/2-1)) ];
else
QM=[eye((M-1)/2) zeros((M-1)/2,1) j*eye((M-1)/2) ; zeros((M-1)/2,1)' sqrt(2) zeros((M-1)/2,1)' ; fliplr(eye((M-1)/2)) zeros((M-1)/2,1) -j*fliplr(eye((M-1)/2)) ];
  QMM1=[eye((M-1)/2) j*eye((M-1)/2) ; fliplr(eye((M-1)/2)) -j*fliplr(eye((M-1)/2))];
end

%Transform array data
Y=kron(QM',QN')*X;

%Compute a basis for the signal subspace
[U,D,S]=svd([real(Y) imag(Y)]);
Es=U(:,1:num_sig);

%Set up more constants
J2M=eye(M);
J2M=J2M(2:M,:);
J4N=eye(N);
J4N=J4N(2:N,:);

K1=real(QMM1'*J2M*QM);
K2=imag(QMM1'*J2M*QM);
Ku1=kron(K1,eye(N));
Ku2=kron(K2,eye(N));
K3=real(QNM1'*J4N*QN);
K4=imag(QNM1'*J4N*QN);
Kv1=kron(eye(M),K3);
Kv2=kron(eye(M),K4);

PSIu=Ku1*Es\Ku2*Es;
PSIv=Kv1*Es\Kv2*Es;

lambda=eig(PSIu+j*PSIv);
u=2*atan(-real(lambda));
v=2*atan(-imag(lambda));

%Compute azimuth and elevation angles of arrival
phi=atan2(v/dy,u/dx);
theta=asin(u./(2*pi*dx*cos(phi)));

%adjust phi if necessary so that 0<=phi<2*pi
for n=1:num_sig
  if phi(n)<0 
   phi(n)=phi(n)+2*pi;
 end
 if phi(n)>=2*pi 
   phi(n)=phi(n)-2*pi;
 end
end

elseif M==1 & N==1
 disp('Array must have more than 1 element')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Equations for Total Least Squares(TLS) ESPRIT
%   See Roy and Kailath's paper(July 1989 IEEE Trans. on
%   Acoustics, SPeech and Sig. Proc.) for a derivation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else  %IF the array is a uniform linear array, then use the TLS ESPRIT

  R=X*X';                 %Compute an estimate of the spatial covariance matrix
  [U,D,V]=svd(R);         %Do eigendecompositon of R

  S=U(:,1:num_sig);       %Compute a basis for the signal subspace

  %Compute the angles of arrival
  theta = S(1:max(N,M)-1,:)\S(2:max(N,M),:);
  w=angle(eig(theta));
  theta=asin(-w/(dx*2*pi));
  phi=[];  %elevation angle is not used in linear array

end


