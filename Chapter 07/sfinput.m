% This m file creat the space factor and save it in sf.m

% Example 1:
% Rectangular SF
% Width: 45 to 135 degrees.
x=0:180;
y=zeros(1,length(x));
for i=1:length(x)
    if (i>=46)&(i<=136)
        y(i)=1;
    end
end
save sf1.m x y -ascii;

% Example 2:
% Triangular SF
% Width: 45 to 135 degrees.
x=0:180;
y=zeros(1,length(x));
for i=1:length(x)
    if (i>=46)&(i<=91)
        y(i)=1/45*x(i)-1;
    elseif (i>=91)&(i<=136)
        y(i)=-1/45*x(i)+3;
    end
end
save sf2.m x y -ascii;