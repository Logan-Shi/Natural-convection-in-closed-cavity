close all;

xcs = linspace(dx/2,Lx-dx/2,nx);
ycs = linspace(dy/2,Ly-dy/2,ny);

pcell = pstar(2:nx+1,2:ny+1);
figure;
[hl,Nul,ql] = calc_heat_transfer(Tcell);
plot(ycs,Nul);
disp(mean(Nul))
title('local Nusselt on left wall');
xlabel('y');
ylabel('Nu');

figure;
quiver(xcs, ycs, ucell', vcell',7);
title('Cell velocity');
xlabel('x');
ylabel('y');
axis equal;
axis([0 Lx 0 Ly]);

figure;
contourf(xcs,ycs,pcell',20);
xlabel('x');
ylabel('y');
axis equal;
axis([0 Lx 0 Ly]);
title("Pressure");

psi = easy_streamfunc(xcs,ycs,nx,ny,dx,dy,ucell,vcell);
figure;
contour(xcs, ycs, psi', 80);
xlabel('x');
ylabel('y');
title('stream-function');
axis equal;
axis([0 Lx 0 Ly]);

figure()
[X,Y,Z] = griddata(xcs,ycs,flipud(rot90(Tcell)),xcs',ycs,'linear');
pcolor(X,Y,Z);
axis equal
axis([0 Lx 0 Ly]);
shading interp;
title('Temperature contour')
colorbar,colormap(jet)