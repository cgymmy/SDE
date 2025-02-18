function tildeu=onestep(u,dt,flux,elementnodes)
    addpath('../SpacialOp');
    u1=u+dt*getRHS(u,flux,elementnodes);
    u2=3/4*u+1/4*u1+1/4*dt*getRHS(u1,flux,elementnodes);
    tildeu=1/3*u+2/3*u2+2/3*dt*getRHS(u2,flux,elementnodes);
    %tildeu=u+dt*getRHS(u,flux,elementnodes);
    return;
end