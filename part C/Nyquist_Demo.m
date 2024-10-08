function Nyquist_Demo
clc;
close all

% Variable
islog = 1;
Yn = 1;
gam = 1;
bet = 0.01;

% Fixed
a1 = 1;
a2 = 1;
m1 = 2;
m2 = 2;
k1 = 2;
k2 = 2*(2*gam+1);
c1 = bet*k1;
c2 = bet*k2;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if islog == 1
    om = 10.^linspace(-1,2,1e4);
else
    om = linspace(0.1,2.5,1e4);
%     om = linspace(0.9,1.15,1e4); % For gam = 0.05
end

Q1_F = a1./(k1 - m1*om.^2 + 1i*c1*om);
Q2_F = a2./(k2 - m2*om.^2 + 1i*c2*om);
if Yn == 1
    Y_F = Q1_F + Q2_F;
else
    Y_F = Q1_F - Q2_F;
end
dQ1_F = 1i*om.*Q1_F;
dQ2_F = 1i*om.*Q2_F;
V_F = 1i*om.*Y_F;

AQ1_F   = 20*log10(abs(Q1_F));
AQ2_F   = 20*log10(abs(Q2_F));
AY_F    = 20*log10(abs(Y_F));

PQ1_F = unwrap(angle(Q1_F));
PQ2_F = unwrap(angle(Q2_F));
PY_F = unwrap(angle(Y_F));

RV_F = real(V_F);
IV_F = imag(V_F);
RdQ1 = real(dQ1_F);
IdQ1 = imag(dQ1_F);
if Yn == 1
    RdQ2 = real(dQ2_F);
    IdQ2 = imag(dQ2_F);
else
    RdQ2 = -real(dQ2_F);
    IdQ2 = -imag(dQ2_F);
end

figure

% Bode plot of amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axA = subplot(2,2,1); hold on; box on;
xlabel('Frequency (rad/sec)')
ylabel('Magnitude (dB)')

YLim = [max(floor(min([AQ1_F,AQ2_F,AY_F])/10)*10,-100), 50];

plot(om,AQ1_F,'r')
plot(om,AQ2_F,'Color',[0,0.65,0])
plot(om,AY_F,'b')

plhAL   = plot([1,1]*om(1),YLim,'k');
plhAQ1  = plot(om(1),AQ1_F(1),'.r','MarkerSize',10);
plhAQ2  = plot(om(1),AQ2_F(1),'.','Color',[0,0.65,0],'MarkerSize',10);
plhAY   = plot(om(1),AY_F(1),'.b','MarkerSize',10);

set(axA,'XLim',om([1,end]), ...
        'YLim',YLim)
if islog == 1
    set(axA,'XScale','log')
end

% Bode plot of phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axP = subplot(2,2,3); hold on; box on;
xlabel('Frequency (rad/sec)')
ylabel('Phase (deg)')

if Yn == 1
    YLim = [-1.1*pi,0.1*pi];
else
    YLim = [-2.1*pi,0.1*pi];
end

plot(om,PQ1_F,'r')
plot(om,PQ2_F,'Color',[0,0.65,0])
plot(om,PY_F,'b')

plhPL   = plot([1,1]*om(1),YLim,'k');
plhPQ1  = plot(om(1),PQ1_F(1),'.r','MarkerSize',10);
plhPQ2  = plot(om(1),PQ2_F(1),'.','Color',[0,0.65,0],'MarkerSize',10);
plhPY   = plot(om(1),PY_F(1),'.b','MarkerSize',10);

set(axP,'XLim',om([1,end]), ...
        'YLim',YLim)
if islog == 1
    set(axP,'XScale','log')
end

% Nyquist plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axN = subplot(1,2,2); hold on; box on; axis square;
xlabel('Real')
ylabel('Imag')

tithN = title(sprintf('omega = %.2f',om(1)));

% plot(RdQ1,IdQ1,'r')
% plot(RdQ2,IdQ2,'Color',[0,0.65,0]);
plot(RV_F,IV_F,'b')

plhNL  = plot(RV_F(1),IV_F(1),'.b','MarkerSize',10);
plhNQ1 = plot([0,RdQ1(1)],[0,IdQ1(1)],'r');
plhNQ2 = plot(RdQ1(1)+[0,RdQ2(1)],IdQ1(1)+[0,IdQ2(1)],'Color',[0,0.65,0]);
plhNY = plot([0,RV_F(1)],[0,IV_F(1)],'b');

% Construct GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UserData = struct(  'axA',      axA, ...
                    'axP',      axP, ...
                    'axN',      axN, ...
                    'Yn',       Yn, ...
                    'gam',      gam, ...
                    'bet',      bet, ...
                    'plhAL',    plhAL, ...
                    'plhAQ1',   plhAQ1, ...
                    'plhAQ2',   plhAQ2, ...
                    'plhAY',    plhAY, ...
                    'plhPL',    plhPL, ...
                    'plhPQ1',   plhPQ1, ...
                    'plhPQ2',   plhPQ2, ...
                    'plhPY',    plhPY, ...
                    'tithN',    tithN, ...
                    'plhNL',    plhNL, ...
                    'plhNQ1',   plhNQ1, ...
                    'plhNQ2',   plhNQ2, ...
                    'plhNY',    plhNY);
set(gcf,'UserData',UserData, ...
        'WindowButtonMotionFcn',    @mouseMove)
end

function mouseMove(object, ~)

% Get mouse position
axA = object.UserData.axA;
XLA = get(axA,'XLim');
YLA = get(axA,'YLim');
CA = get(axA, 'CurrentPoint');
xA = CA(1,1);
yA = CA(1,2);

axP = object.UserData.axP;
XLP = get(axP,'XLim');
YLP = get(axP,'YLim');
CP = get(axP, 'CurrentPoint');
xP = CP(1,1);
yP = CP(1,2);

% If mouse inside of axes
if xA > XLA(1) && xA < XLA(2) && yA > YLA(1) && yA < YLA(2)
    om = xA;
elseif xP > XLP(1) && xP < XLP(2) && yP > YLP(1) && yP < YLP(2)
    om = xP;
else
    return
end

% Variable
Yn = object.UserData.Yn;
gam = object.UserData.gam;
bet = object.UserData.bet;

% Fixed
a1 = 1;
a2 = 1;
m1 = 2;
m2 = 2;
k1 = 2;
k2 = 2*(2*gam+1);
c1 = bet*k1;
c2 = bet*k2;

Q1_F = a1./(k1 - m1*om.^2 + 1i*c1*om);
Q2_F = a2./(k2 - m2*om.^2 + 1i*c2*om);
if Yn == 1
    Y_F = Q1_F + Q2_F;
else
    Y_F = Q1_F - Q2_F;
end
dQ1_F = 1i*om.*Q1_F;
dQ2_F = 1i*om.*Q2_F;
V_F = 1i*om.*Y_F;

AQ1_F   = 20*log10(abs(Q1_F));
AQ2_F   = 20*log10(abs(Q2_F));
AY_F    = 20*log10(abs(Y_F));

PQ1_F = unwrap(angle(Q1_F));
PQ2_F = unwrap(angle(Q2_F));
PY_F = unwrap(angle(Y_F));

RV_F = real(V_F);
IV_F = imag(V_F);
RdQ1 = real(dQ1_F);
IdQ1 = imag(dQ1_F);
if Yn == 1
    RdQ2 = real(dQ2_F);
    IdQ2 = imag(dQ2_F);
else
    RdQ2 = -real(dQ2_F);
    IdQ2 = -imag(dQ2_F);
end

% Bode plot of amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(object.UserData.plhAL,'XData',[1,1]*om)
set(object.UserData.plhAQ1,'XData',om,'YData',AQ1_F)
set(object.UserData.plhAQ2,'XData',om,'YData',AQ2_F)
set(object.UserData.plhAY,'XData',om,'YData',AY_F)

% Bode plot of phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(object.UserData.plhPL,'XData',[1,1]*om)
set(object.UserData.plhPQ1,'XData',om,'YData',PQ1_F)
set(object.UserData.plhPQ2,'XData',om,'YData',PQ2_F)
set(object.UserData.plhPY,'XData',om,'YData',PY_F)

% Nyquist plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(object.UserData.tithN,'String',sprintf('omega = %.2f',om))
set(object.UserData.plhNL,'XData',RV_F,'YData',IV_F)
set(object.UserData.plhNQ1,'XData',[0,RdQ1],'YData',[0,IdQ1])
set(object.UserData.plhNQ2,'XData',RdQ1+[0,RdQ2],'YData',IdQ1+[0,IdQ2])
set(object.UserData.plhNY,'XData',[0,RV_F],'YData',[0,IV_F])
end