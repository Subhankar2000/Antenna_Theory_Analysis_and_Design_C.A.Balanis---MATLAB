function y = Subroutine(x)

global G a b Lambda f;
    F1=((sqrt(2.*x)-b/Lambda).^2).*(2.*x-1);
    F2=(((G/(2*pi))*sqrt(3/(2*pi)).*sqrt(1./x)-a/Lambda).^2);
    F3=(((G^2)/(6*pi^3))./x-1);
    y=F1-F2.*F3;
