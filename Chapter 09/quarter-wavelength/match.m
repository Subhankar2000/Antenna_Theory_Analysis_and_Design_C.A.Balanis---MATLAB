function result=match(Z0,Rl,Xl);

% ESTIMATION OF THE CLOSEST POINT TO THE LOAD (s WAVELENGTHS)   
%
% THAT THE NEW LOAD IS ONLY REAL AND DERIVATION OF THIS NEW LOAD (new Rl)
%
% INPUT PARAMETERS ARE THE TRANSMISSION LINE CHARACTERISTIC IMPEDANCE (Z0)
%
% AND THE IMPEDANCE LOAD, REAL (Rl) AND IMAGINARY PART (Xl)

Zl=Rl+j*Xl;

polynom=[Z0*Xl (-Z0^2+(abs(Zl))^2) -Z0*Xl];

roots_polynom=roots(polynom);

dist=(1/(2*pi))*atan(roots_polynom);

accept_dist=find(dist>=0);

new_dist=min(dist(accept_dist));

tankl=tan(2*pi*new_dist);

new_Zl=Z0*(Zl+j*Z0*tankl)/(Z0+j*Zl*tankl);

new_Rl=real(new_Zl);

result=[new_dist new_Rl];
