function [theta,theta_samp,b,varargout]=wood(str,d,len,Nel,sfaf_mode,filename,pstring,Nsamp,Ntheta,theta,SFAF);

close all;

if ~isempty(pstring),
   px=inline(vectorize(pstring));
else
   px=[];
end;

if (mod(Nsamp,2)==0),   % even number    
    m=1:(Nsamp/2);
    theta_samp=acos([(2*m-1)/(2*len) (1-2*m)/(2*len)])*180/pi;
    costhm=[(2*m-1)/(2*len) (1-2*m)/(2*len)];

%    costhm(find(abs(imag(theta_samp))>=1e-2))=[];
%    costhm(abs(costhm)>1)=[];

    theta_samp(abs(imag(theta_samp))>=1e-2)=[];
    theta_samp=real(theta_samp);

    costhm=cos(theta_samp*pi/180);

    switch sfaf_mode,
       case 1  % file 
          b=interp1(theta,SFAF,theta_samp,'linear');  % linear interpolation
       case 2  % function 
          b=px(theta_samp*pi/180);
    end;

else     % odd number
    m=-(Nsamp-1)/2:(Nsamp-1)/2;
    theta_samp=acos(m/len)*180/pi;
%     theta_samp(abs(imag(theta_samp))>=1e-2)=[];
%     costhm=m/len;
% %    costhm(find(abs(imag(theta_samp))>=1e-2))=[];
% 
     theta_samp=real(theta_samp);
    costhm=cos(theta_samp*pi/180);

    switch sfaf_mode,
       case 1  % file 
          b=interp1(theta,SFAF,theta_samp,'linear');  % linear interpolation
       case 2  % function 
          b=px(theta_samp*pi/180);
    end;
  
end;  

SF_rec=zeros(size(theta));
AF_rec=zeros(size(theta));


for m=1:length(costhm),
   if str=='Space Factor',
      SF_rec=SF_rec+b(m)*sin(pi*len*(cos(theta*pi/180)-costhm(m)))./ ...
             (pi*len*(cos(theta*pi/180)-costhm(m)));
   else
      AF_rec=AF_rec+b(m)*sin(Nel*pi*d*(cos(theta*pi/180)-costhm(m)))./ ...
             (Nel*sin(pi*d*(cos(theta*pi/180)-costhm(m))));
   end;
end;

if str=='Space Factor',
   zs=linspace(-len/2,len/2,length(theta));
   Is=zeros(size(zs));

   for m=1:length(theta_samp),
      Is=Is+1/len*b(m)*exp(-j*2*pi*zs*costhm(m));
   end;

   varargout={zs,Is,SF_rec};

else

   if mod(Nel,2)==0,   % even samples
      za=((1:Nel/2)-0.5)*d;
      za=[-fliplr(za) za];
   else   % odd samples
      za=(-(Nel-1)/2:(Nel-1)/2)*d;
   end;

   Ia=zeros(size(za));
   
   for m=1:Nel,
      Ia(m)=Ia(m)+1/Nel*sum(b.*exp(-j*2*pi*za(m)*costhm));
   end;

   varargout={za,Ia,AF_rec};

end;

