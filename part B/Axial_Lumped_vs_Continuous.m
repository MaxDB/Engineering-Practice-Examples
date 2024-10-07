clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
N = 10;     % Number of elements

rho = 8000; % Density [kg/m^3]
A = 1e-4;   % Area [m^2]
E = 200e9;  % Young's modulus [Pa]
l = 1;      % Length [m]

% Compute natural frequencies and modeshapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

om_L = zeros(2,N);
om_C = zeros(2,N);
for n = 1:N
    L = l/n;    % Length of element

    % Lumped mass local matrices
    M_L_local = rho*A*L/2*[ 1,  0;
                            0,  1];
    K_L_local = E*A/L*[	1,  -1;
                        -1, 1];

    % Continous mass local matrices
    M_C_local = rho*A*L/6*[ 2,  1;
                            1,  2];
    K_C_local = E*A/L*[ 1,  -1;
                        -1, 1];

    % Global matrices
    M_L = zeros(n+1,n+1);
    K_L = zeros(n+1,n+1);
    M_C = zeros(n+1,n+1);
    K_C = zeros(n+1,n+1);

    for i = 1:n
        idx = [i,i+1];
        M_L(idx,idx) = M_L(idx,idx) + M_L_local;
        K_L(idx,idx) = K_L(idx,idx) + K_L_local;
        M_C(idx,idx) = M_C(idx,idx) + M_C_local;
        K_C(idx,idx) = K_C(idx,idx) + K_C_local;
    end

    % Find eigenvalues/eigenvectors
    [Vec_L, Val_L] = eig(M_L\K_L);
    [Vec_C, Val_C] = eig(M_C\K_C);

    % Rearrange modes (in ascending order)
    Val_L = diag(Val_L);
    [Val_L,idx] = sort(Val_L);
    Vec_L = Vec_L(:,idx);

    Val_C = diag(Val_C);
    [Val_C,idx] = sort(Val_C);
    Vec_C = Vec_C(:,idx);

    % Scale modeshapes
    Vec_L = bsxfun(@rdivide,Vec_L,max(abs(Vec_L)));
    Vec_C = bsxfun(@rdivide,Vec_C,max(abs(Vec_C)));
    
    % Record natural frequencies
    if n == 1
        om_L(1,n) = sqrt(Val_L(2));
        om_C(1,n) = sqrt(Val_C(2));
    else
        om_L(:,n) = sqrt(Val_L(2:3));
        om_C(:,n) = sqrt(Val_C(2:3));
    end
end

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True natural frequencies
om2_T = pi/l*sqrt(E/rho);
om3_T = 2*pi/l*sqrt(E/rho);

figure

% Plot of second natural frequencies
om2_L = om_L(1,:);
om2_C = om_C(1,:);

subplot(2,2,1)
hold on; box on;
xlim([1,N])
plot([1,N],[1,1]*om2_T,'k')
plot(1:N,om2_L,'.-')
plot(1:N,om2_C,'.-')
xlabel('Number of elements')
ylabel('\omega_2')
legend( 'True','Lumped','Continuous', ...
        'Location','SE')

% Plot of third natural frequencies
om3_L = om_L(2,2:end);
om3_C = om_C(2,2:end);

subplot(2,2,2)
hold on; box on;
xlim([2,N])
plot([2,N],[1,1]*om3_T,'k')
plot(2:N,om3_L,'.-')
plot(2:N,om3_C,'.-')
xlabel('Number of elements')
ylabel('\omega_3')
legend( 'True','Lumped','Continuous', ...
        'Location','SE')

% Plot of errors of second natural frequencies
om2_L_ERR = log10(abs(om2_L - om2_T)/om2_T);
om2_C_ERR = log10(abs(om2_C - om2_T)/om2_T);

subplot(2,2,3)
hold on; box on;
xlim([1,N])
plot(1:N,om2_L_ERR,'.-')
plot(1:N,om2_C_ERR,'.-')
xlabel('Number of elements')
ylabel('log_{10} of error of \omega_2')
legend( 'Lumped','Continuous', ...
        'Location','NE')

% Plot of errors of second natural frequencies
om3_L_ERR = log10(abs(om3_L - om3_T)/om2_T);
om3_C_ERR = log10(abs(om3_C - om3_T)/om2_T);

subplot(2,2,4)
hold on; box on;
xlim([2,N])
plot(2:N,om3_L_ERR,'.-')
plot(2:N,om3_C_ERR,'.-')
xlabel('Number of elements')
ylabel('log_{10} of error of \omega_3')
legend( 'Lumped','Continuous', ...
        'Location','NE')