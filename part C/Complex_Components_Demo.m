function Complex_Components_Demo

close all

t = linspace(0,1,1e3);
om = 2*pi*2;
A = 1;
phi = 0;

y = A*cos(om*t + phi);

yp = A/2*exp(+1i*(om*t + phi));
ym = A/2*exp(-1i*(om*t + phi));

subplot(1,2,1)
hold on; box on;
plot(t,y,'k')
plh0 = plot(t(1),y(1),'.k','MarkerSize',10);
xlabel('t')
ylabel('y')

subplot(1,2,2)
hold on; box on;
plot3(t,real(yp),imag(yp),'b')
plot3(t,real(ym),imag(ym),'r')
plhp = plot3(t(1),real(yp(1)),imag(yp(1)),'.b','MarkerSize',10);
plhm = plot3(t(1),real(ym(1)),imag(ym(1)),'.r','MarkerSize',10);
xlabel('t')
ylabel('real(y)')
zlabel('imag(y)')
legh = legend('1/2e^{+i(\omega_2 t + \phi)}','1/2e^{-i(\omega_2 t + \phi)}');
set(legh,'Position',get(legh,'Position')+[0.075,0.075,0,0])
view([-50,15])

disp('Press any key to continue...')
pause

SlowMo = 8;
tic
while toc < t(end)*SlowMo
    [~,idx] = min(abs(t*SlowMo - toc));
    set(plh0,'XData',t(idx),'YData',y(idx))
    set(plhp,'XData',t(idx),'YData',real(yp(idx)),'ZData',imag(yp(idx)))
    set(plhm,'XData',t(idx),'YData',real(ym(idx)),'ZData',imag(ym(idx)))
    drawnow
end