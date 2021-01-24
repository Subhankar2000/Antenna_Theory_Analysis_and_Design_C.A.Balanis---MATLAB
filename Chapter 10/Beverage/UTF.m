function y = UTF(TH)

global l K

y = (sin(TH).^2).*( sin( pi*l.*(cos(TH)-K)) ).^2./(cos(TH)-K).^2;
