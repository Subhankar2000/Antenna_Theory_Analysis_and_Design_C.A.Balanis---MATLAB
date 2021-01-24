function [result]=fact(x);

for i=1:length(x),
   result(i)=factorial(x(i))*(x(i)>0)+(x(i)==0);
end;   






