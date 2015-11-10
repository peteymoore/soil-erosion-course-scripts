% Green-Ampt excess rainfall model
% P. Moore - Sept. 2015
%   based on Examples 3.2 amd 3.3 from Julien, 2002
clear all; %close all;

% ====== Soil properties =========

p0 = 0.501;         % -, total porosity
p0e = 0.486;        % -, drainable porosity
Se = 0.3;           % -, initial effective saturation <==== ANTECEDENT MOISTURE
psi = 16.68;        % cm, suction head
K = 0.65;           % cm/hr, saturated hydraulic conductivity
hstor = 0.75;       % cm, depression storage capacity
dtheta0 = p0e*(1-Se);   % -, moisture capacity

% ====== Rainfall time series =========

t = 0:15:180;       % min, time elapsed
thr = t/60;         % hr, converted
idt = [0 0.25 0.3 0.43 0.66 1.55 4.85 0.91 0.51 0.36 0.28 0.23 0.2];
                    % cm, rainfall during time increment
nt = length(t);
dt = thr(2)-thr(1);
                    
% ======= Pre-allocate arrays ================

fp = zeros(nt,1); dFa = zeros(nt,1); Fp = zeros(nt,1); dFp = zeros(nt,1);
Sstor = zeros(nt,1); exc = zeros(nt,1); runoff = zeros(nt,1);
dSstor = zeros(nt,1); dexc = zeros(nt,1); idt_plus_Sstor = zeros(nt,1);

% * dFa is actual cumulative infiltration. fp is the potential infiltration 
% rate and Fp and dFp the potential cumulative infiltration and step change 
% (independent of rainfall intensity). All are functions of time, so 
% require arrays of length nt. Matlab (and other languages) run more
% efficiently if variable sizes don't change within loops.

% Sstor is surface detention storage in micro-depressions, up to hstor. exc
% is the excess that can become runoff once depression storage is full.

% ======= Set up initial conditions =====================
%    * first two rows are ICs; all pre-filled with zeros

Fp(1) = 1e-5; % need nonzero value since Fp is in denominator of fp eqn.


% ======= Solution loop =================================

for i = 2:nt
    fp(i) = K*(1 + (psi*dtheta0)/Fp(i-1));      % pot. inf. rate  <======
    dFp(i) = fp(i)*dt;                          % pot. cumul. inf.
    idt_plus_Sstor(i) = idt(i)+Sstor(i-1);      % water needing entry
    dFa(i) = min(dFp(i),idt_plus_Sstor(i));     % actual inf. increment
    Fp(i) = Fp(i-1) + dFa(i);                   % actual infilt. <========
    dexc(i) = idt_plus_Sstor(i) - dFa(i);       % excess rainfall
    Sstor(i) = min(dexc(i),hstor);              % depression storage
    runoff(i) = dexc(i) - Sstor(i);             % actual runoff
end

% ======= Plot resulting infiltration & runoff ==========

%figure;
subplot(2,2,1); stairs(t,idt,'b-'); ylabel('rainfall, cm'); hold on
subplot(2,2,2); plot(t,dFa*4,'r-o'); ylabel('infiltration rate, cm/hr'); hold on
subplot(2,2,3); plot(t,Fp,'b-o'); ylabel('infiltration, cm'); xlabel('time, min'); hold on
subplot(2,2,4); stairs(t,runoff,'r-'); ylabel('runoff, cm'); xlabel('time, min'); hold on