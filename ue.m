function uv=ue(z)
x=z(1); y=z(2);
%uv=exp(x+y);                       Neumann +Dirichlet
%uv=x*y;                            Neumann + Dirichlet
%uv=x*(x-1)*(y-1)*y;                only Dirichlet

% uv=x*y+x+y;
if (x==1)||(x==0)
    uv=pi/2;
end
if(y==0)||(y==1)
    uv=0;
end