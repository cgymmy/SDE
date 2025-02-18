%Basis Function on given intervals!
%deg:the degree of function
%dir:the direvative of function
%[left,right]:the support interval
function y=basis1D(x,deg,dir,left,right)
    if deg==0
        if dir==0
            y=1;
        elseif dir==1
            y=0;
        else
            output("Error derivative");
        end
    elseif deg==1
        if dir==0
            y=1/(right-left)*(x-(left+right)/2);
        elseif dir==1
            y=1/(right-left);
        else
            output("Error derivative");
        end
    elseif deg==2
        h=right-left;
        if dir==0
            y=1/2*(3*((2*x-left-right)/h).^2-1);
        elseif dir==1
            y=6*(2*x-left-right)/(h.^2);
        else
            output("Error derivative")
        end
    end
    return;
end