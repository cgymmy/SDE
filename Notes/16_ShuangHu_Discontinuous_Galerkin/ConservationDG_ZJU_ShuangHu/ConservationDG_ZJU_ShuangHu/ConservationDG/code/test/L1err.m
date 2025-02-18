function e=L1err(res,truesol,elementnodes,T)
    n=size(elementnodes,2)-1;
    e=0;
    for i=1:n
        func=@(x)abs(assembleu(res,x,elementnodes(i),elementnodes(i+1),i)-truesol(x,T));
        e=e+myintegral(func,elementnodes(i),elementnodes(i+1));
    end
    return;
end