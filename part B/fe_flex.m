
%fixed/fixed with support at midspan.

n=6 %number of elements NEEDS TO BE EVEN

%M=rho AL/420 Mstar with y defined in terms of y_1, theta_1*L etc:
Mstarlocal=[156 22 54 -13;22 4 13 -3;54 13 156 -22; -13 -3 -22 4];
%K=EI/L Kstar with y defined in terms of y_1, theta_1*L etc:
Kstarlocal=[12 6 -12 6;6 4 -6 2;-12 -6 12 -6;6 2 -6 4];

Mstar=zeros(2*n+2,2*n+2);
Kstar=zeros(2*n+2,2*n+2);
for i=1:n %each element in turn
    Mstar((2*i-1):(2*i+2),(2*i-1):(2*i+2))=Mstar((2*i-1):(2*i+2),(2*i-1):(2*i+2))...
        +Mstarlocal;
    Kstar((2*i-1):(2*i+2),(2*i-1):(2*i+2))=Kstar((2*i-1):(2*i+2),(2*i-1):(2*i+2))...
        +Kstarlocal;
end  

%apply boundary conditions

Mstar=Mstar([2:n n+2:end-2 end],[2:n n+2:end-2 end]);
Kstar=Kstar([2:n n+2:end-2 end],[2:n n+2:end-2 end]);

%eigenvalues/vectors: eig(X) solves (X-lambda I)U=0
% normalised eigenvectors (in terms of nodes of U) are the columns in w 
% and corresponding to eigenvalues in leading diagonal of lambda
[w,lambda]=eig(inv(Mstar)*Kstar);

%omega=omegastar*sqrt(EI/(rho A))/l^2  noting L=l/n where l is the total beam length
% therefore omegastar=sqrt(lambda*420)*n^2 since omega^2=420/L^4 EI/(rhoA)lambda^2
omegastar=sqrt(diag(lambda)*420)*n^2;
%sort into order of modes
[omegastar,shuffle]=sort(omegastar);
w=w(:,shuffle);

%mode shapes
%add zeros to w corresponding to the fixed points:
wl=size(w,1);
w=[zeros(1,wl);w(1:(wl-1)/2,:);zeros(1,wl);w((wl+1)/2:wl-1,:);zeros(1,wl);w(wl,:)];
for mode=[1:3]
    figure(mode)
    hold off
    xoverL=0:.01:1;
        yglobal=zeros(1,100*n+1);
for element=1:n
        y=(1-3*xoverL.^2+2*xoverL.^3)*w(2*element-1,mode)+...
        (xoverL-2*xoverL.^2+xoverL.^3)*w(2*element,mode)+...
        (3*xoverL.^2-2*xoverL.^3)*w(2*element+1,mode)+...
        (-xoverL.^2+xoverL.^3)*w(2*element+2,mode);
    yglobal(((element-1)*100+1):(element*100))=y(1:100);
end    
    %plot((element-1)/n:.01/n:(element)/n,y/max(abs(y)))
    plot((0:0.01/n:1),yglobal/max(abs(yglobal)))
    hold on
    plot((0:1/n:1),yglobal(1:100:end)/max(abs(yglobal)),'o')
end
plot(0:.001:1,-sin(4*pi*(0:.001:1)),'g')
figure(2)
figure(1)
plot(0:.001:1,sin(2*pi*(0:.001:1)),'g')
omegastar(1:3)
    