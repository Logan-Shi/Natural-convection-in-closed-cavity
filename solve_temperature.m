% ## Copyright (C) 2013 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## calc_temperature

% ## Author: logan-shi <loganshi@sjtu.edu.cn>
% ## Created: 2020-06-07

function T = solve_temperature(u,v,T)

simple_globals;

relax = Trelax;

SIZE = (nx-1) * (ny-1);
A = spalloc(SIZE,SIZE,5*SIZE);
% A = zeros(SIZE,SIZE);

% vectorized version
dte = alpha_heat * dy / dx;
dtw = alpha_heat * dy / dx;
dtn = alpha_heat * dx / dy;
dts = alpha_heat * dx / dy;

fte = dy / 2 * (u(2:nx,1:ny-1)+u(2:nx,2:ny));
ftw = dy / 2 * (u(1:nx-1,1:ny-1)+u(1:nx-1,2:ny));
ftn = dx / 2 * (v(1:nx-1,2:ny)+v(2:nx,2:ny));
fts = dx / 2 * (v(1:nx-1,1:ny-1)+v(2:nx,1:ny-1));
%Ò»½×Ó­·çÊ½
ate = dte + max(-fte,0);
atw = dtw + max(ftw,0);
atn = dtn + max(-ftn,0);
ats = dts + max(fts,0);
atp = ate+atw+atn+ats+fte-ftw+ftn-fts+rho*dx*dy/dt;

% RHS
% rhs = rho*dx*dy/dt * T(2:nx,2:ny)...
%     +(1-relax)/relax * atp .* T(2:nx,2:ny);
rhs = zeros(nx-1,ny-1);
% % E
% rhs(nx,:) = rhs(nx,:) + ate(nx,:) .* T(nx+1,2:ny+1);
% ate(nx,:) = 0;
% % W
% rhs(1,:) = rhs(1,:) + atw(1,:) .* T(1,2:ny+1);
% atw(1,:) = 0;
% % N
% rhs(:,ny) = rhs(:,ny) + atn(:,ny) .* T(2:nx+1,ny);
% atn(:,ny) = 0;
% % S
% rhs(:,1) = rhs(:,1) + ats(:,1) .* T(2:nx+1,1);
% ats(:,1) = 0;

idx = 0;
for j = 2:ny
    for i = 2:nx
        idx = idx + 1;

        S = ats(i-1,j-1)*T(i,j-1);
        W = atw(i-1,j-1)*T(i-1,j);
        E = ate(i-1,j-1)*T(i+1,j);
        N = atn(i-1,j-1)*T(i,j+1);
        rhs(i-1,j-1) = S+W+E+N;
        A(idx,idx) = atp(i-1,j-1);
    end
end

rhs = reshape(rhs,[SIZE,1]);

% solve
sol = A \ rhs;
Tmp = reshape(sol,nx-1,ny-1);
T(2:nx,2:ny) = Tmp;%flipud(rot90(Tmp));

% u = blkdiag(0,u,0);
% v = blkdiag(0,v,0);
% TT = T;
% for i = 2:nx+1
%     for j = 2:ny+1
%         T(i,j) = TT(i,j)+dt*((-0.25)*(((u(i,j)+u(i+1,j))*(TT(i,j)+...
%             TT(i+1,j))-(u(i-1,j)+u(i,j))*TT(i-1,j)+TT(i,j))/dx+...
%             ((v(i,j)+v(i,j+1))*(TT(i,j)+TT(i,j+1))-(v(i,j)+v(i,j-1))*(TT(i,j)+TT(i,j-1)))/dy)...
%             +(alpha_heat/dx/dy)*(TT(i+1,j)+TT(i-1,j)+TT(i,j+1)+TT(i,j-1)-4*TT(i,j)));
%     end
% end

% T = T(2:nx+1,2:ny+1);
% A = zeros(nx*ny);
% b = zeros(nx*ny,1);
% for i = 2:nx-1
%     for j = 2:ny-1
%         A(i+(j-1)*nx,i-1+(j-1)*nx) = 1/dy/dy;
%         A(i+(j-1)*nx,i+1+(j-1)*nx) = 1/dx/dx;
%         A(i+(j-1)*nx,i+(j-1)*nx) = -2/dy/dy-2/dx/dx;
%         A(i+(j-1)*nx,i+(j-2)*nx) = 1/dx/dx;
%         A(i+(j-1)*nx,i+j*nx) = 1/dx/dx;
%         dTdx = (T(i-1,j)+T(i+1,j)-2*T(i,j))/dx/dx;
%         dTdy = (T(i,j-1)+T(i,j+1)-2*T(i,j))/dy/dy;
%         b(i+(j-1)*nx) = (u(i,j)*dTdx+v(i,j)*dTdy)/alpha_heat;
%     end
% end
% 
% for i = 1:nx
%     A(i,i) = 1;
%     b(i) = T(i,1);
%     A(i+(ny-1)*nx,i+(ny-1)*nx) = 1;
%     b(i+(ny-1)*nx) = T(i,ny);
% end
% 
% for j = 2:ny-1
%     A(1+(j-1)*nx,1+(j-1)*nx) = 1;
%     b(1+(j-1)*nx) = T(1,j);
%     A(j*nx,j*nx) = 1;
%     b(j*nx) = T(nx,j);
% end
% sol = A\b;
% T = reshape(sol,nx,ny);
% T = blkdiag(0,T,0);
% enforce BC
T = apply_theta_bc(T);

return
end