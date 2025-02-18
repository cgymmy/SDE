function tildeu=postprocess(u,elementnodes)
    addpath('../SpacialOp')
    m=size(u,1);
    n=size(u,2);
    tildeu=zeros(m,n);
    if m==3
        for i=1:n
            % get values
            [x3,x1]=limiter(u,i,elementnodes);
            x2=assembleu(u,(elementnodes(i)+elementnodes(i+1))/2,elementnodes(i),elementnodes(i+1),i);
            % interpolation
            tildeu(1,i)=(x1+x3)/6+2/3*x2;
            tildeu(2,i)=x3-x1;
            tildeu(3,i)=(x1+x3)/3-2/3*x2;
        end
    elseif m==1
        error("Don't need postprocess!");
    elseif m==2
        for i=1:n
            % get values
            [x2,x1]=limiter(u,i,elementnodes);
            % interpolation
            tildeu(1,i)=(x1+x2)/2;
            tildeu(2,i)=x2-x1;
        end        
    end
    return;
end