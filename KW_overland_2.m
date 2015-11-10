% KW_overland_2.m - overland kinematic flow solver
% -------------------------------------------------------------------------
%       P. Moore, 9/2015
%
% A 1D model for overland flow with the kinematic wave approximation.
% Solution uses the MacCormack finite difference scheme as described in
% Huang and Lee, 2013, Journal of Hydrology 489(238-245).
% 
% Can be adapted to arbitrary slope morphology, roughness distribution, and
% storm properties.
% -------------------------------------------------------------------------
close all; clear all;

% ========= generate hillslope geometry ==========
L = 60;                         % m, length
dx = 3;                         % m, spatial step
x = 0:dx:L;
nx = length(x);
theta = 0.1;                    % m/m, slope
z0 = 50;                        % m, summit elevation
dz = dx*theta;
z = z0:-dz:(z0-(nx-1)*dz);

% ======== flow resistance/rating equation =================
mann = 0.1;                     % Manning's n
% below is an anonymous function, called in the for loop that follows
disch = @(H) 1/mann*theta^(1/2)*H^(5/3);

% ======== set up time stepping ===============================
dt = 2;
t = 0:dt:45*60;
nt = length(t);

% ========= input excess rainfall (default or output from Green-Ampt_1.m ==
i = 2/100/(60*60);         % excess rainfall, input as cm/hr [m/s]
I = i.*ones(nt,nx);
I(500:end,:) = 0;

% ========= pre-allocate arrays ============================
q = zeros(nt,nx); h = zeros(nt,nx); dhadj = zeros(1,nx);
dh = zeros(1,nx);

% ========= upslope boundary condition, initial condition ==============
h(1,:) = 0;
h(:,1) = 0;

% ======= main solution loop ==========================================
%
% here, dhadj is an adjustment to the increment in flow depth that helps to
% stabilize the solution. The magnitude of dhadj depends on lambda, which
% is a measure of the celerity of a kinematic wave compared with the
% ratio of discretization parameters dx/dt

for j = 2:nt
    for k = 2:nx-1
        dh(j,k) = dt/dx*(disch(h(j-1,k-1)) - disch(h(j-1,k))) +...
            0.5*(I(j,k)+I(j-1,k))*dt;
        lam(k) = max(0,sqrt(h(j-1,k)*9.81)-(dx/dt));
        dhadj(k) = (dh(j,k)+lam(k)*(dt/dx)*dhadj(k-1))/(1+lam(k)*(dt/dx));
        h(j,k) = h(j-1,k) + dhadj(k);
        q(j,k) = disch(h(j,k));
    end
    h(j,end) = h(j,end-1);        % forces downslope b.c. dh/dx = dq/dx = 0
end

% ========= Vlsualize the result! ========================= 
subplot(2,1,1);
plot(x,h(1:20:end,:),'b-');
xlabel('distance downslope, m'); ylabel('mean flow depth, m');

subplot(2,1,2);
plot(t(1:20:end)./60,q(1:20:end,3),'b-o'); hold on
plot(t(1:20:end)./60,q(1:20:end,ceil(nx/2)),'m-o');
plot(t(1:20:end)./60,q(1:20:end,nx-1),'r-o');
xlabel('time,min'); ylabel('unit discharge, m^2/s');
legend('upslope','midslope','footslope');