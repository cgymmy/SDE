%ODE:U_{t}=\mathcal{L}(U)
%This function we assemble the RHS vector U.
%input:the initial function u, the flux function, the element nodes
%dummy:the right boundary as dummy node
%output:\mathcal{L}(u)
function res=getRHS(u,flux,elementnodes)
    m=size(u,1);
    n=size(u,2);
    res=zeros(m,n);
    dummy=elementnodes(n+1);
    for element=1:n
        for j=1:m
            left=elementnodes(element);
            right=elementnodes(element+1);
            if element==1
                newelem=n;
            else
                newelem=element-1;
            end
            res(j,element)=InnerProduct(flux,u,j-1,left,right,element)...
                -LaxFriedrichesFlux(flux,u,element,elementnodes)*basis1D(right,j-1,0,left,right)...
                +LaxFriedrichesFlux(flux,u,newelem,elementnodes)*basis1D(left,j-1,0,left,right);
            massfunc=@(x)basis1D(x,j-1,0,left,right).^2;
            mass=myintegral(massfunc,left,right);
            res(j,element)=res(j,element)/mass;
        end
    end
    return;
end