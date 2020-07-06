load('20Ra1e5.mat')
load('20Ra1e5.mat')
Nul20 = Nul;
ycs20 = ycs;
load('40Ra1e5.mat')
load('40Ra1e5.mat')
Nul40 = Nul;
ycs40 = ycs;
load('60Ra1e5Nul.mat')
load('60Ra1e5y.mat')
Nul60 = Nul;
ycs60 = ycs;
load('70Ra1e5Nul.mat')
load('70Ra1e5y.mat')
Nul70 = Nul;
ycs70 = ycs;
load('75Ra1e5.mat')
load('75Ra1e5.mat')
Nul75 = Nul;
ycs75 = ycs;
load('80Ra1e5Nul.mat')
load('80Ra1e5y.mat')
Nul80 = Nul;
ycs80 = ycs;
load('85Ra1e5.mat')
load('85Ra1e5.mat')
Nul85 = Nul;
ycs85 = ycs;
figure;
hold on
plot(ycs20,Nul20);
plot(ycs40,Nul40);
plot(ycs60,Nul60);
plot(ycs70,Nul70);
plot(ycs75,Nul75);
plot(ycs80,Nul80);
plot(ycs85,Nul85);
legend('20','40','60','70','75','80','85')
title('local Nusselt on left wall');
xlabel('y');
grid on
ylabel('Nu');