function y = UT(TH,PHI)

global h l K

y =  ( sin( 2*pi*h.*sin(TH).*cos(PHI) ).^2 ).* ... 
     (sin(TH).^2).*( sin( pi*l.*(cos(TH)-K)) ).^2./(cos(TH)-K).^2;
