function init=getinitval(elementnodes,order)
    Nx=size(elementnodes,2)-1;
    if order==0
        init=initval(elementnodes(1:Nx));
    elseif order==1
        init=zeros(2,Nx);
        for j=1:Nx
            left=initval(elementnodes(j));
            right=initval(elementnodes(j+1));
            init(1,j)=(left+right)/2;
            init(2,j)=right-left;
        end
    elseif order==2
        interpnodes=zeros(1,2*Nx+1);
        for i=1:2:2*Nx+1
            interpnodes(i)=elementnodes((i+1)/2);
        end
        for i=2:2:2*Nx
            interpnodes(i)=(interpnodes(i+1)+interpnodes(i-1))/2;
        end
        interpval=initval(interpnodes(1:2*Nx+1));
        init=zeros(3,Nx);
        for i=1:Nx
            x1=interpval(2*i-1);
            x2=interpval(2*i);
            x3=interpval(2*i+1);
            init(1,i)=(x1+x3)/6+2/3*x2;
            init(2,i)=x3-x1;
            init(3,i)=(x1+x3)/3-2/3*x2;
        end
    end
end