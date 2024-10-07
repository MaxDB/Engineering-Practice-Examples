function Signal_From_Fourier_Demo

close all

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

% Plot of signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axX = subplot(3,1,1); hold on; box on;

YLimX = [-2,+2];
t = linspace(0,1,1e3);
x = zeros(size(t));
plhX = plot(t,x);
xlim([0,1])
ylim(YLimX)
xlabel('time')
ylabel('x')

% Cos amp axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axC = subplot(3,1,2); hold on; box on;

Ac = zeros(1,10);
XLimC = [-0.25,length(Ac)-0.75];
YLimC = [-1,1];

plot(XLimC,[0,0],'k')
plhC1 = plot(0:length(Ac)-1,Ac,'.b','MarkerSize',16);
plhC2 = plot(0:length(Ac)-1,Ac,'or','MarkerSize',4);

% xlabel('Frequency')
ylabel('Cos amplitde')
xlim(XLimC)
ylim(YLimC)
set(gca,'XTickLabel',[])

% Sin amp axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axS = subplot(3,1,3); hold on; box on;

As = Ac;
XLimS = XLimC;
YLimS = YLimC;

plot(XLimS,[0,0],'k')
plhS1 = plot(0:length(As)-1,As,'.b','MarkerSize',16);
plhS2 = plot(0:length(As)-1,As,'or','MarkerSize',4);

xlabel('Frequency')
ylabel('Sin amplitde')
xlim(XLimS)
ylim(YLimS)
set(gca,'XTickLabel',sprintf('%i\\omega\n',0:9))

% Construct GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UserData = struct(  'plhX',     plhX, ...
                    'plhC',     [plhC1,plhC2], ...
                    'plhS',     [plhS1,plhS2], ...
                    'axX',      axX, ...
                    'axC',      axC, ...
                    'axS',      axS, ...
                    'YLimX',    YLimX, ...
                    'XLimC',    XLimC, ...
                    'YLimC',    YLimC, ...
                    'XLimS',    XLimS, ...
                    'YLimS',    YLimS, ...
                    't',        t);
set(gcf,'UserData',UserData, ...
        'WindowButtonMotionFcn',    @mouseMove, ...
        'WindowButtonDownFcn',      @mouseClick)
end

function mouseMove(object, ~)

% Get mouse position
axC = object.UserData.axC;
CC = get(axC, 'CurrentPoint');
xC = CC(1,1);
yC = CC(1,2);

axS = object.UserData.axS;
CS = get(axS, 'CurrentPoint');
xS = CS(1,1);
yS = CS(1,2);

plhC = object.UserData.plhC;    % Get Fourier plot handles
Ac = get(plhC(1),'YData');      % Get fixed Fourier values

plhS = object.UserData.plhS;    % Get Fourier plot handles
As = get(plhS(1),'YData');      % Get fixed Fourier values

% Get axis limits
XLimC = object.UserData.XLimC;
YLimC = object.UserData.YLimC;
XLimS = object.UserData.XLimS;
YLimS = object.UserData.YLimS;

% If mouse inside of cos axis
if xC > XLimC(1) && xC < XLimC(2) && yC > YLimC(1) && yC < YLimC(2)
    k = round(xC);   % Get Fourier frequency for mouse position
    
    % Set new temp. Fourier value
    Ac(k+1) = yC;
    set(plhC(2),'YData',Ac)
else
    set(plhC(2),'YData',Ac) % Set circle positions to dot position
end

% If mouse inside of sin axis
if xS > XLimS(1) && xS < XLimS(2) && yS > YLimS(1) && yS < YLimS(2)
    k = round(xS);   % Get Fourier frequency for mouse position
    
    % Set new temp. Fourier value
    As(k+1) = yS;
    set(plhS(2),'YData',As)
else
    set(plhS(2),'YData',As) % Set circle positions to dot position
end

% Recompute and set signal
plhX = object.UserData.plhX;
t = object.UserData.t;
x = zeros(size(t));
for i = 1:10
    x = x + Ac(i)*cos((i-1)*2*pi*t) + ...
            As(i)*sin((i-1)*2*pi*t);
end
set(plhX,'YData',x)

YLimX = object.UserData.YLimX;
MMx = 1.1*[min(x),max(x)];

if MMx(1) < YLimX(1)
    YLimX(1) = MMx(1);
end

if MMx(2) > YLimX(2)
    YLimX(2) = MMx(2);
end

axX = object.UserData.axX;
set(axX,'ylim',YLimX)
end

function mouseClick(object, ~)


% Get mouse position
axC = object.UserData.axC;
CC = get(axC, 'CurrentPoint');
xC = CC(1,1);
yC = CC(1,2);

axS = object.UserData.axS;
CS = get(axS, 'CurrentPoint');
xS = CS(1,1);
yS = CS(1,2);

plhC = object.UserData.plhC;    % Get Fourier plot handles
Ac = get(plhC(1),'YData');      % Get fixed Fourier values

plhS = object.UserData.plhS;    % Get Fourier plot handles
As = get(plhS(1),'YData');      % Get fixed Fourier values

% Get axis limits
XLimC = object.UserData.XLimC;
YLimC = object.UserData.YLimC;
XLimS = object.UserData.XLimS;
YLimS = object.UserData.YLimS;

% If mouse inside of cos axis
if xC > XLimC(1) && xC < XLimC(2) && yC > YLimC(1) && yC < YLimC(2)
    k = round(xC);   % Get Fourier frequency for mouse position
    
    % Set new temp. Fourier value
    Ac(k+1) = yC;
    set(plhC(1),'YData',Ac)
end

% If mouse inside of sin axis
if xS > XLimS(1) && xS < XLimS(2) && yS > YLimS(1) && yS < YLimS(2)
    k = round(xS);   % Get Fourier frequency for mouse position
    
    % Set new temp. Fourier value
    As(k+1) = yS;
    set(plhS(1),'YData',As)
end

% Recompute and set signal
plhX = object.UserData.plhX;
t = object.UserData.t;
x = zeros(size(t));
for i = 1:10
    x = x + Ac(i)*cos((i-1)*2*pi*t) + ...
            As(i)*sin((i-1)*2*pi*t);
end
set(plhX,'YData',x)

YLimX = object.UserData.YLimX;
MMx = 1.1*[min(x),max(x)];

if MMx(1) < YLimX(1)
    YLimX(1) = MMx(1);
end

if MMx(2) > YLimX(2)
    YLimX(2) = MMx(2);
end

axX = object.UserData.axX;
set(axX,'ylim',YLimX)
end
