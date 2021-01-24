function [prec,HPBW_app,tm]=find_hpbw(theta,U)

% U is radiation intensity vector
% *******************************
U=U(theta<=pi);
theta=theta(theta<=pi);

[Um,im]=max(U);
th_3db=[];
prec=1e-3;

while isempty(th_3db),
   th_3db=theta(abs(U-0.5*Um)<=prec);
%   th_low=max(th_3db(th_3db<=theta(im)));
%   th_high=min(th_3db(th_3db>=theta(im)));
   if isempty(th_3db)|length(th_3db)<2,
      prec=prec*10;
   end;     
end;   

if theta(im)<min(th_3db)|theta(im)>max(th_3db),
   HPBW_app=(min(th_3db)+abs(pi-max(th_3db)))*180/pi;      
else
   th_low=max(th_3db(th_3db<theta(im)));
   th_high=min(th_3db(th_3db>theta(im)));
   HPBW_app=abs(th_low-th_high)*180/pi;   
end;

tm=theta(im);
