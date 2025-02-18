%function:InnerProduct
%Input:function u,flux function,the degree of test function testdeg, the
%jth element:elementid
%interval [left,right], the jth element
%Output:\int_{left}^{right}f(u)v_x\dif x.
function res=InnerProduct(flux,u,testdeg,left,right,elementid)
    ux=@(x)assembleu(u,x,left,right,elementid);
    func=@(x)flux(ux(x))*basis1D(x,testdeg,1,left,right);
    %res=integral(func,left,right,'ArrayValued',true);
    res=myintegral(func,left,right);
end