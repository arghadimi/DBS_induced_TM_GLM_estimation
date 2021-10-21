function I=SteadyStateCurrent(param,F_dbs)
f=param.f;
U=param.U;
F=param.F;
D=param.D;
t_syn=param.t_syn;
u=(f+U*(1-f)*(1-exp(-1./(F_dbs*F))))/(1-(1-f)*exp(-1./(F_dbs*F))); % TM descrete
% u=(U+(f-U)*exp(-1./(F_dbs*F)))/(1-(1-f)*exp(-1./(F_dbs*F)));% Costa et al.
r=((1-exp(-1./(F_dbs*D)))/(1-(1-u)*exp(-1./(F_dbs*D))));%(u*(1-f)+f)*
% Tau_syn=3;

I=r*u/(1-exp(-1./(F_dbs*t_syn)));
