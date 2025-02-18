%Newton's method
%solve the nonlinear function
%to get the exact solution of Burgers equation
function y=exactBurgers(x,t)
    if x<=pi+0.5*t
        y=1;
    else
        y=-1;
    end
    tol=1e-10;
    while(abs(y-sin(x-y*t-0.5*t))>tol)
        y=y-(y-sin(x-y*t-0.5*t))/(1+t*cos(x-y*t-0.5*t));
    end
    y=y+0.5;
    return;
end