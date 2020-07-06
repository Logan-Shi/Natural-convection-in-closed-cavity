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

% ## SimpleBox

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-02

clc;
% clf;
clear all;

simple_globals;

% initialization
Lx = 1;
Ly = 1;
% ncell = 40;
ncell = 70;
nx = ncell;
ny = ncell;
dx = Lx / nx;
dy = Ly / ny;

relax = 0.5;
urelax = relax;
vrelax = relax;
prelax = 0.8;
Trelax = 0.8;
accuracy = 4e-3;

ULid = 0;
rho = 1.165;   % density,kg/m^3
visc = 18.6e-6;   % viscosity,kg/m/s
kvis = 16e-6; % kinematic viscosity,m^2/s
T0 = 300; %K
T_H = 320;
T_L = 280;
Gry = 9.81e-6/3/1.19; %m/s^2
alpha_heat = 22.9e-6;     % heat transmission coeffient,m^2/s
lambda = 0.0267; % heat conduction,W/m/K

% dt = 0.025 * Re * dx^2;
% dt = 0.001;
dt = Inf;

% storage
umac = zeros(nx+1,ny+2);
vmac = zeros(nx+2,ny+1);
ustar = zeros(nx+1,ny+2);
vstar = zeros(nx+2,ny+1);
p = zeros(nx+2,ny+2);
pstar = zeros(nx+2,ny+2);
pdash = zeros(nx+2,ny+2);
T = zeros(nx+1,ny+1)+T0;
T = apply_theta_bc(T);

umac = apply_umac_bc(umac);
vmac = apply_vmac_bc(vmac);

uold = umac;
vold = vmac;
pold = p;

AUp = zeros(nx+1,ny+2);
AVp = zeros(nx+2,ny+1);

% computation begins
Pr = kvis/alpha_heat;
Gr = Gry*1/T0*(T_H-T_L)*Ly^3/kvis/kvis;
Ra = Gr*Pr;
disp(['closed cavity, natural convection, Gr=', num2str(Gr)]);
disp(['closed cavity, natural convection, Ra=', num2str(Ra)]);
disp(['grid size: ' num2str(ncell) '*' num2str(ncell)])
debug = 0;
max_steps = 5000;
hist_acc = 0;
counter = 0;
for istep = 1:max_steps
    close all;
    Told = T;
    uold = umac;
    vold = vmac;
    pold = pstar;
    
    ustar = solve_ustar(umac,vmac,uold,vold,pstar);
    vstar = solve_vstar(umac,vmac,uold,vold,pstar,T,T0);
    
    [umac,vmac,pstar] = solve_pressure(ustar,vstar,pstar);
    
    ucell = 0.5 * (umac(1:nx,2:ny+1) + umac(2:nx+1,2:ny+1));
    vcell = 0.5 * (vmac(2:nx+1,1:ny) + vmac(2:nx+1,2:ny+1));


    [T] = solve_temperature(ucell,vcell,T);
    Tcell = (T(1:nx,1:ny)+T(1:nx,2:ny+1)+T(2:nx+1,1:ny)+T(2:nx+1,2:ny+1))/4;
    % compute intermeidate velocity
    
    if debug
    figure()
    xcs = linspace(dx/2,Lx-dx/2,nx);
    ycs = linspace(dy/2,Ly-dy/2,ny);
    quiver(xcs, ycs, ucell', vcell',3);
    title('Cell velocity');
    xlabel('x');
    ylabel('y');
    axis equal;
    axis([0 Lx 0 Ly]);
    
    figure()
    [X,Y,Z] = griddata(xcs,ycs,flipud(rot90(Tcell)),xcs',ycs,'v4');
    pcolor(X,Y,Z);
    axis equal
    axis([0 Lx 0 Ly]);
    shading interp;
    title('Temperature contour')
    colorbar,colormap(jet)
    end

    if (mod(istep,10) == 0)
        
        ucorr = umac - uold;
        vcorr = vmac - vold;
        Tcorr = T - Told;
        % pcorr = pstar - pold;
        tol_abs = accuracy*norm(umac,inf);
        corr = norm(ucorr,inf);
        counter = counter+1;
        hist_acc(counter) = corr/tol_abs*accuracy;
        disp(['step=', int2str(istep), ...
            ', corr=', num2str(hist_acc(counter))]);
        
        if (corr < tol_abs)
            break;
        end
    end
end

close all;

xcs = linspace(dx/2,Lx-dx/2,nx);
ycs = linspace(dy/2,Ly-dy/2,ny);

pcell = pstar(2:nx+1,2:ny+1);
figure;
[hl,Nul,ql] = calc_heat_transfer(Tcell);
plot(ycs,Nul);
title('local heat convection coefficient on left wall');
xlabel('y');
ylabel('hl');

figure;
quiver(xcs, ycs, ucell', vcell',3);
title('Cell velocity');
xlabel('x');
ylabel('y');
axis equal;
axis([0 Lx 0 Ly]);

% figure;
% contourf(xcs,ycs,pcell',20);
% xlabel('x');
% ylabel('y');
% axis equal;
% axis([0 Lx 0 Ly]);
% title("Pressure");
% 
% psi = easy_streamfunc(xcs,ycs,nx,ny,dx,dy,ucell,vcell);
% figure;
% contour(xcs, ycs, psi', 80);
% xlabel('x');
% ylabel('y');
% title('stream-function');
% axis equal;
% axis([0 Lx 0 Ly]);

figure()
[X,Y,Z] = griddata(xcs,ycs,flipud(rot90(Tcell)),xcs',ycs,'linear');
pcolor(X,Y,Z);
axis equal
axis([0 Lx 0 Ly]);
shading interp;
title('Temperature contour')
colorbar,colormap(jet)




