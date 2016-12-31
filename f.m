function fv=f(z);
x=z(1); y=z(2);
%fv=-exp(x+y);                                 Dirichlet and Neumann
%%%%%%%%fv=(x^2-x)*(y^2-y)-2*(x^2+y^2-x-y);    Homogenous dirichlet
%fv=0;
fv=0;