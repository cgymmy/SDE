function [adjuminus,adjuplus]=limiter(u,eid,elementnodes)
    %% Manage the element id.
    n=size(elementnodes,2)-1;
    if eid+1>n
        righteid=1;
    else
        righteid=eid+1;
    end
    if eid>1
        lefteid=eid-1;
    else
        lefteid=n;
    end
    %% Assemble functions and derive ubar
    func=@(x)assembleu(u,x,elementnodes(eid),elementnodes(eid+1),eid);
    rightfunc=@(x)assembleu(u,x,elementnodes(righteid),elementnodes(righteid+1),righteid);
    leftfunc=@(x)assembleu(u,x,elementnodes(lefteid),elementnodes(lefteid+1),lefteid);
    bar=myintegral(func,elementnodes(eid),elementnodes(eid+1))/(elementnodes(eid+1)-elementnodes(eid));
    rightbar=myintegral(rightfunc,elementnodes(righteid),elementnodes(righteid+1))/(elementnodes(righteid+1)-elementnodes(righteid));
    leftbar=myintegral(leftfunc,elementnodes(lefteid),elementnodes(lefteid+1))/(elementnodes(lefteid+1)-elementnodes(lefteid));
    %% Get tilde and double tilde
    tilde=func(elementnodes(eid+1))-bar;
    doubletilde=bar-func(elementnodes(eid));
    rightdelta=rightbar-bar;
    leftdelta=bar-leftbar;
    %% get the minmod function
    if sign(tilde)==sign(leftdelta) && sign(leftdelta)==sign(rightdelta)
        tildemod=min([abs(leftdelta),abs(rightdelta),abs(tilde)])*sign(tilde);
    else
        tildemod=0;
    end
    if sign(doubletilde)==sign(leftdelta) && sign(leftdelta)==sign(rightdelta)
        doubletildemod=min([abs(leftdelta),abs(rightdelta),abs(doubletilde)])*sign(doubletilde);
    else
        doubletildemod=0;
    end
    %% get the result
    adjuminus=bar+tildemod;
    adjuplus=bar-doubletildemod;
end