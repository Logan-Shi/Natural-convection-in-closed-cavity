function [hl,Nul,ql] = calc_heat_transfer(T)
simple_globals;
dT = T(1,:)-T_H;
ql = -lambda*dT/(dx/2);
hl = ql/(T_H-T_L);
Nul = hl*Ly/lambda;
end
% function [hl,Nul,ql] = calc_heat_transfer(T)
% simple_globals;
% dT = T(floor(1+size(T,1)/2),:)-T(floor(1+size(T,1)/2)-1,:);
% ql = -lambda*dT/dx;
% hl = ql/(T_H-T_L);
% Nul = hl*Ly/lambda;
% end