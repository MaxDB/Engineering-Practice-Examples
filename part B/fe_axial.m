% axial beam vibration
close all
clc
%number of elements
n=20;
%modes to be plotted
m=[1:4];
%end conditions: 0=free 1=fixed
lefthandside=1;
righthandside=1;

%local matricies
%Mlocal=(rho AL/6)Mstarlocal
Mstarlocal=[2 1;1 2];
%Klocal=(EA/L)Kstarlocal
Kstarlocal=[1 -1;-1 1];

%global
Mstar=zeros(n+1,n+1);
Kstar=zeros(n+1,n+1);
for i=1:n %each element in turn
    Mstar(i:i+1,i:i+1)=Mstar(i:i+1,i:i+1)+Mstarlocal;
    Kstar(i:i+1,i:i+1)=Kstar(i:i+1,i:i+1)+Kstarlocal;
end  
%end conditions
if lefthandside==1
    Mstar=Mstar(2:end,2:end);    
    Kstar=Kstar(2:end,2:end);
end
if righthandside==1
    Mstar=Mstar(1:end-1,1:end-1);    
    Kstar=Kstar(1:end-1,1:end-1);
end

% (-omega^2 M+K)U=0 gives (inv(Mstar)Kstar-(rho L^2/(6E))omega^2 I)U=0
% or (X-lambda I)U=0 with X=inv(Mstar)Kstar and lambda=(rho L^2/(6E))omega^2

%eigenvalues/vectors: eig(X) solves (X-lambda I)U=0
% normalised eigenvectors (in terms of nodes of U) are the columns in w 
% and corresponding to eigenvalues in leading diagonal of lambda
[w,lambda]=eig(inv(Mstar)*Kstar);
%omega=omegastar*sqrt(E/rho)/l  noting L=l/n where l is the total beam length
% therefore omegastar=sqrt(rho/E)*L*n*sqrt(lambda*6E/(rho L^2))=sqrt(lambda*6)*n
omegastar=sqrt(diag(lambda)*6)*n;
%sort into order of modes
[omegastar,shuffle]=sort(omegastar);
w=w(:,shuffle);

%note if free-free first mode appears as rigid body motion - remove:
if [lefthandside==0 & righthandside==0]
    omegastar=omegastar(2:end);
    w=w(:,2:end);
end


%add in zero nodes to eigenvector equation
if lefthandside==1
    w=[zeros(1,size(w,2));w];    
end
if righthandside==1
    w=[w;zeros(1,size(w,2))];    
end

%plotting
no_of_modes=n+1-lefthandside-righthandside-[lefthandside==0 & righthandside==0];
a=find(m>=no_of_modes);
m(a)=no_of_modes;

x=0:(1/n/10):1;
if [lefthandside==1 & righthandside==1]
    for i=m
        subplot(size(m,2),1,find(m==i))
        node_pts=11:10:(size(x,2)-1); %locations in x vector of non-zero node points
        plot([0 x(node_pts) 1],w(:,i),'.-')
        hold on
        shape=sin(pi*i*x);
        normamp=sqrt(sum(shape(node_pts).^2)); %normalised at same points as approx sol.
        signamp=sign(shape(node_pts(1)))*sign(w(2,i)); %ensures signs match
        plot(x,signamp*shape/normamp,'g') %exact solution
    end
    omegastarexact=pi*(1:no_of_modes)'; %exact solution
    omegastar_approx_exact=[omegastar omegastarexact]
elseif [lefthandside==0 & righthandside==0]
    for i=m
        subplot(size(m,2),1,find(m==i))
        node_pts=1:10:size(x,2); %locations in x vector of non-zero node points
        plot(x(node_pts),w(:,i),'.-')
        hold on
        shape=-cos(pi*i*x);
        normamp=sqrt(sum(shape(node_pts).^2)); %normalised at same points as approx sol.
        signamp=sign(shape(node_pts(1)))*sign(w(1,i)); %ensures signs match
        plot(x,signamp*shape/normamp,'g') %exact solution
    end
    omegastarexact=pi*(1:no_of_modes)'; %exact solution
    omegastar_approx_exact=[omegastar omegastarexact]
else
    plot((1:1:n+1)',w(:,1:m),'.-')
    xlabel('node')
    ylabel('normalised modeshape')
    omegastar
end


      
        
        

