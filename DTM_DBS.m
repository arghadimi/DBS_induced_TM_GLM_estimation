function S=DTM_DBS(param,F_DBS,pulse_num)
    f=param.f;
    U=param.U;
    F=param.F;
    D=param.D;
    t_syn=param.t_syn;
    
    u=U*(1-f)+f;
    R=1;
    I_sparse=u*R;

    dT=1/F_DBS;
    for i=1:pulse_num-1
        [u,R,I_sparse(i+1)]=DBS_DTM(f,U,F,D,t_syn,u,R,I_sparse(i),dT);
    end
    S=I_sparse;
end
function [u_new,R_new,I]=DBS_DTM(f,U,F,D,syn,u_old,R_old,I,dT)

    R_new=1+(R_old*(1-u_old)-1)*exp(-dT/D);

    u_new=U+f*(1-U)+(1-f)*(u_old-U)*exp(-dT/F); %% TM decret
    % u_new=U+(u_old+f*(1-u_old)-U)R*exp(-dT/F); %% Costa  et al.
    I=I*exp(-dT/syn)+R_new*u_new;
end