function Basis=RCosBasis(hist)
dt=hist.dt;

% fie=hist.fie;
RCos= @ (a,c,fie,t) .5*cos(max(fie-pi,min(fie+pi,a'*log(dt*t+c)))-fie)+.5;
t=1:hist.delay;



Basis=[];
for c=hist.c
    for a=hist.a
        for fie=hist.fie;
            Basis(end+1,:)=RCos(a,c,fie,t);
        end
    end
end