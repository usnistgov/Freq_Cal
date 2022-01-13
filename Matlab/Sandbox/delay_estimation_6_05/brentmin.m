function [xmin, fx, iter]=brentmin(ax,bx,cx,f,tol,itmax)
%[xmin, fx, iter]=brentmin(ax,bx,cx,f,tol,ITMAX)
%This solves for the minimum of function f. The minimum xmin must be 
%bracketed by ax and cx, and bx must be a point on this interval with a
%value less than at either end-point. The algorithm ends when a relative
%tolerance tol is reached or the maximum number of iterations itmax is
%reached.  This optionally returns the function value at the minimum, fx, 
%and the number of iterations, iter.
%
%This code is based on the "brent" function in Press's Numerical Recipes in C.
%
% Copyright Travis Wiens t.wiens@usask.ca 2015

if nargin<5||isempty(tol)
    tol=1e-5;
end

if nargin<6||isempty(itmax)
    itmax=100;
end

CGOLD=0.3819660;
ZEPS=1.0e-10;


e=0.0;

if ax < cx
    a=ax;
else
    a=cx;
end
if ax>cx
    b=ax;
else
    b=cx;
end

x=bx;
w=bx;
v=bx;


fx=f(x);
fw=fx;
fv=fx;

d=nan;

for iter=1:itmax
    xm=0.5*(a+b);
    tol1=tol*abs(x)+ZEPS;
    tol2=2.0*tol1;
    if (abs(x-xm) <= (tol2-0.5*(b-a)))
        xmin=x;
        return
    end
    if (abs(e) > tol1)
        r=(x-w)*(fx-fv);
        q=(x-v)*(fx-fw);
        p=(x-v)*q-(x-w)*r;
        q=2.0*(q-r);
        if (q > 0.0)
            p = -p;
        end
        q=abs(q);
        etemp=e;
        e=d;
        if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            if (x >= xm)
                e=a-x;
            else
                e=b-x;
            end
            d=CGOLD*e;
        else
            d=p/q;
            u=x+d;
            if (u-a < tol2 || b-u < tol2)
                d=tol1*sign(xm-x);
            end
        end
        
    else
        if (x >= xm)
            e=a-x;
        else
            e=b-x;
        end
        d=CGOLD*e;
    end
    if (abs(d) >= tol1)
        u=x+d;
    else
        u=x+tol1*sign(d);
    end
    fu=f(u);
    if (fu <= fx)
        if (u >= x)
            a=x;
        else
            b=x;
        end
        %SHFT(v,w,x,u)
        v=w;
        w=x;
        x=u;
        %SHFT(fv,fw,fx,fu)
        fv=fw;
        fw=fx;
        fx=fu;
    else
        if (u < x)
            a=u;
        else
            b=u;
        end
        if (fu <= fw || w == x)
            v=w;
            w=u;
            fv=fw;
            fw=fu;
        elseif (fu <= fv || v == x || v == w)
            v=u;
            fv=fu;
        end
    end
end
xmin=x;
end
