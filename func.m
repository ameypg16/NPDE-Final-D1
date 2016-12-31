function f=func(z)
x=z(1); y=z(2);
d=0.06;
if (y==0)||(y==1)
    if (x<=d)&&(x>=0)
        f=x/d;
    elseif (x<=1-d)&&(x>=d)
        f=1;
    elseif (x>=1-d)&&(x<=1)
        f=(1-x)/d;
    end
elseif (x==0)||(x==1)
    if (y<=d)&&(y>=0)
        f=y/d;
    elseif (y<=1-d)&&(y>=d)
        f=1;
    elseif (y>=1-d)&&(y<=1)
        f=(1-y)/d;
    end
else
    f=1;
end
f=abs(f);

    
    