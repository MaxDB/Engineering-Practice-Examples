function Modeshape_Demo

close all

m1 = 5; % Mass 1
m2 = 5; % Mass 2

k1 = 50; % Spring 1
k2 = 50; % Spring 2
k3 = 10; % Spring 3

T = 10; % Time of simulation

FPS = 30;   % Target FPS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create mass and stiffness matrices
M = [   m1, 0;
        0,  m2];
K = [   k1+k2,  -k2;
        -k2,    k2+k3];

% Find the Lambda and Phi matrices
[Phi,Lam] = eig(M\K);

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

blue_col = [0.3,0.3,1];

% Masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
hold on
axis equal off
set(gca,'XLim',[-1,+1]*2.1)

% Create masses
rh1 = rectangle('position',[-1.5,-0.5,1,1],...
                'curvature',[0.05,0.05],...
                'FaceColor',blue_col);
rh2 = rectangle('position',[+0.5,-0.5,1,1],...
                'curvature',[0.05,0.05],...
                'FaceColor',blue_col);

% x1 vs x2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
hold on
box on
axis equal square
set(gca,'XLim',[-1,+1]*0.5, ...
        'YLim',[-1,+1]*0.5, ...
        'XTick',[], ...
        'YTick',[])

plot([0,0],[-0.5,+0.5],'k')
plot([-0.5,+0.5],[0,0],'k')

modeH1 = plot([-1,+1]*Phi(1,1),[-1,+1]*Phi(2,1),'--k','Visible','off');
modeH2 = plot([-1,+1]*Phi(1,2),[-1,+1]*Phi(2,2),'--k','Visible','off');

doth = plot(0,0,'.b','MarkerSize',15);
xlabel('x_1')
ylabel('x_2')

% Add checkbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axPos = get(gca,'Position');
uicontrol(  'Style','PushButton', ...
            'Units','Normalized', ...
            'Position',[axPos(1),0.9,0.2,0.05], ...
            'String','Modeshapes On', ...
            'Callback', @ButtonPush, ...
            'UserData', [modeH1,modeH2])

% Construct GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UserData = struct(  'handles',  [doth,rh1,rh2], ...
                    'Phi',      Phi, ...
                    'Lam',      Lam, ...
                    'T',        T, ...
                    'FPS',      FPS);
set(gcf,'UserData',UserData, ...
        'WindowButtonMotionFcn',    @mouseMove, ...
        'WindowButtonDownFcn',      @mouseClick)
end

function mouseMove(object, ~)

UserData = object.UserData;
doth = UserData.handles(1);
rh1 = UserData.handles(2);
rh2 = UserData.handles(3);

% Get click position
C = get(gca, 'CurrentPoint');

x1 = C(1,1);
x2 = C(1,2);

if abs(x1) > 0.5
    x1 = sign(x1)*0.5;
end
if abs(x2) > 0.5
    x2 = sign(x2)*0.5;
end

set(doth,'XData',x1,'YData',x2)
set(rh1,'position',[-1.5 + x1,-0.5,1,1])
set(rh2,'position',[+0.5 + x2,-0.5,1,1])
end

function mouseClick(object, ~)

% Switch off mouse motion
set(gcf,'WindowButtonMotionFcn', [])
    
UserData = object.UserData;
doth = UserData.handles(1);
rh1 = UserData.handles(2);
rh2 = UserData.handles(3);

Phi = UserData.Phi;
Lam = UserData.Lam;
T   = UserData.T;
FPS = UserData.FPS;

% Get click position
C = get(gca, 'CurrentPoint');

x1 = C(1,1);
x2 = C(1,2);

if abs(x1) > 0.5
    x1 = sign(x1)*0.5;
end
if abs(x2) > 0.5
    x2 = sign(x2)*0.5;
end

Q0 = Phi\[x1;x2];

t = linspace(0,T,T*FPS);

% Find the linear natural frequencies
om1 = sqrt(Lam(1,1));
om2 = sqrt(Lam(2,2));

% Time-domain displacements
q = [   Q0(1)*cos(om1*t);
        Q0(2)*cos(om2*t)];
x = Phi*q;

plot(x(1,1),x(2,1),'.b','MarkerSize',10)
plh = plot(x(1,1),x(2,1),'-');

% Find if modeshapes are switched on
FigKids = get(gcf,'Children');
MS_on = strcmp(get(FigKids(1),'String'),'Modeshapes Off');

if MS_on == 1
    xq1 = Phi(:,1)*q(1,:);
    xq2 = Phi(:,2)*q(2,:);
    
    pdh1 = plot(xq1(1,1),xq1(2,1),'.r','MarkerSize',8);
    pdh2 = plot(xq2(1,1),xq2(2,1),'.r','MarkerSize',8);
    plh1 = plot([xq1(1,1),x(1,1),xq2(1,1)],[xq1(2,1),x(2,1),xq2(2,1)],'-r');
end

% Animate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
while toc < t(end)
    [~,i] = min(abs(t - toc));
    
    set(plh,'XData',x(1,1:i),'YData',x(2,1:i))
    set(doth,'XData',x(1,i),'YData',x(2,i))
    set(rh1,'position',[-1.5 + x(1,i),-0.5,1,1])
    set(rh2,'position',[+0.5 + x(2,i),-0.5,1,1])
    
    if MS_on == 1
        set(pdh1,   'XData',xq1(1,i), ...
                    'YData',xq1(2,i))
        set(pdh2,   'XData',xq2(1,i), ...
                    'YData',xq2(2,i))
        set(plh1,   'XData',[xq1(1,i),x(1,i),xq2(1,i)], ...
                    'YData',[xq1(2,i),x(2,i),xq2(2,i)])
    end
    
    drawnow
end
delete(pdh1)
delete(pdh2)
delete(plh1)

% Switch mouse motion on again
set(gcf,'WindowButtonMotionFcn', @mouseMove)
end

function ButtonPush(source,~)
modeH1 = source.UserData(1);
modeH2 = source.UserData(2);

if strcmp(source.String,'Modeshapes On')
    set(modeH1,'Visible','On')
    set(modeH2,'Visible','On')
    set(source,'String','Modeshapes Off')
else
    set(modeH1,'Visible','Off')
    set(modeH2,'Visible','Off')
    set(source,'String','Modeshapes On')
end
end