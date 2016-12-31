function [uxv,uyv]=uxe(z)
x=z(1); y=z(2);
%uxv=exp(x+y);
%uyv=exp(x+y);
%uxv=(2*x-1)*(y^2-y);
%uyv=(2*y-1)*(x^2-x);
%uxv=0;uyv=0;
uxv=y+1;
uyv=x+1;