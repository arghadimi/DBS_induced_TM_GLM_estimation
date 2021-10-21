function [p]=ScatterPlot(I_bin, dt, up, down, color)
    tspk=find(I_bin>0)*dt;
    p_tot=plot([tspk;tspk],[down*ones(size(tspk));up*ones(size(tspk))],color);
    p=p_tot(1);
end
    