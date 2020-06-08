function [hl,Nul,ql] = calc_heat_transfer(T)
simple_globals;
dT = diff(T);
ql = -lambda*dT(2,:)/dx;
ql = ql(2:end-1);
hl = ql/(T_H-T_L);
Nul = hl*Ly/lambda;
end