function Nyquist_Freq_Demo
clc; close all

om = 10;
f = 10;

t = linspace(0,2*pi,1+100*f);
x = cos(om*t);

N = ceil(2*om/f-1);
y = cos((om-N*f)*t);

figure
hold on; box on;
plot(t([1,end]),[0,0],'Color',[1,1,1]*0.5)
for i = 1:f-1
    plot([1,1]*t(100*i+1),[-1,+1]*1.1,'Color',[1,1,1]*0.5)
end

plot(t,x,'b','Linewidth',1)
for i = 0:f
    plot(t(100*i+1),x(100*i+1),'ob')
end
plot(t,y,'--r','Linewidth',1)

xlim(t([1,end]))
ylim([-1,+1]*1.1)