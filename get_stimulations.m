function I_stim=get_stimulations(param)
    mode=param.mode;
    dt=param.dt;
    T=param.T;
    if strcmpi(mode,'DBS')
        F_DBS=param.F_DBS;

        t_stim=[0];
        I_stim=zeros(1,round(T/dt));%I_rand;%filter(1,[1 -0.9],I_rand);%
        for F_dbs=F_DBS
            T_dbs=1/(dt*F_dbs);
            t_stim=floor((t_stim(end)+1/dt):T_dbs:(t_stim(end)+(1+T)/dt));
            I_stim(t_stim)=DTM_DBS(param.TM,F_dbs,length(t_stim))/DTM_DBS(param.TM,10,1);
        end
    elseif strcmpi(mode,'rand')
        I_rand = randn(1,round(T/dt));
        I_stim=filter(1,[1 -0.9],I_rand);
    elseif strcmpi(mode,'Poisson') || strcmpi(mode,'Poiss')
        I_stim=HomoPoiss(T,dt,param.F_poiss);
    else 
        error('The mode is not valid!')
    end
    
end

function spike=HomoPoiss(T,dt,rr)
    % just for simulation of a spike train
    N = round(T/dt);
    spike=zeros(1,N);
    t=linspace(0,T,N);

    for i=6:N
        if rand<1-exp(-rr*dt) && sum(spike(i-5:i-1))==0;
           spike(i)=1;
        end
    end
end
    