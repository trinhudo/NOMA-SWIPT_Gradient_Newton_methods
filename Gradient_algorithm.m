close all
clear all
%% Simulation parameters
%
PS_dB   = 50;               % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF     = 10;               % S-F distance
dSN     = 3;
dNF     = dSF - dSN;
L       = 1e3;              % path-loss at reference distance
%
lSN     = L*dSN^-3;         % lambda
lSF     = L*dSF^-3;
lNF     = L*dNF^-3;
%
eta     = 0.7;              % energy conversion coefficient
pN      = 0.1;              % power allocation coefficient
pF      = 1 - pN;
RthN    = 1;                % target data rate of User N bits/s/Hz
RthF    = 1;                % target data rate of User N bits/s/Hz
g1 = 2^(2*RthN) - 1;        % gamma_2 for User F
g2 = 2^(2*RthF) - 1;        % gamma_2 for User F
%%
for ss = 1:length(PS_dB)
    % GRADIENT parameters
    k = 0; % step index
    rho = 0.5; % initial value of rho
    delta_k = 1; % initial value of d_k
    stop_th = 10^-5; % stopping threshold
    step_k = 10^-1; % step
    %
    a1 = (1-rho)*pF*PS(ss)/(naN + ncN); % omitting (1-rho) at noise power
    a2 = (1-rho)*pN*PS(ss)/(naN + ncN);
    b1 = pF * PS(ss) / (naF + ncF);
    b2 = pN * PS(ss) / (naF + ncF);
    c  = eta*rho*PS(ss)/(naF + ncF);
    mu_a = g2/(a1-a2*g2);
    mu_b = g2/(b1-b2*g2);
    theta = pF/pN;
    % for optimize
    r1 = -RthN;
    r2 = -RthF;
    snrN = PS(ss)/(naN + ncN);
    snrF = PS(ss)/(naF + ncF);
    nu1 = - g1/(pF-pN*g2)/snrN/lSN;
    nu2 = -g1/pN/snrN/lSN;
    ka = - g2/(pF-pN*g2)/snrN/lSN;
    kb = - mu_b/lSF;
    kc = - g2/lSN/lNF/eta/snrF;
    zeta = 1 - exp(-mu_b/lSF);
    save_rho = rho;
    save_f = [];
    %
    G_k = 1; % temporary G_k
    % SEARCHING OPTIMAL VALUE PROCESS
%     while delta_k > stop_th % old (wrong) stopping criterion
    while norm(G_k) > stop_th % update 01-Aug-2017, use L2-norm
        % Grandient of f(rho)
        if g2/(pF-pN*g2) >= g1/pN % omitting the condition of g2<theta
            f_rho = r1*exp(nu1/(1-rho)) ...
                + r2*exp(ka/(1-rho)+kb) ...
                + r2*zeta*exp(ka/(1-rho)) ...
                + r2*zeta*kc/rho*igamma(0,-ka/(1-rho));
            %             save_f = [save_f, f_rho];
            %
            G_k = r1*nu1*exp(nu1/(1-rho)) / ((1-rho)^2) ...
                + r2*ka*exp(ka/(1-rho)+kb) / ((1-rho)^2) ...
                + r2*ka*zeta*exp(ka/(1-rho)) / ((1-rho)^2) ...
                - r2*kc*zeta*exp(ka/(1-rho))/(1-rho)/rho ...
                - r2*kc*zeta/(rho^2)*igamma(0,-ka/(1-rho));
        else
            f_rho = r1*exp(nu2/(1-rho)) ...
                + r2*exp(ka/(1-rho)+kb) ...
                + r2*zeta*exp(ka/(1-rho)) ...
                + r2*zeta*kc/rho*igamma(0,-ka/(1-rho));
            %             save_f = [save_f, f_rho];
            %
            G_k = r1*nu2*exp(nu2/(1-rho)) / ((1-rho)^2) ...
                + r2*ka*exp(ka/(1-rho)+kb) / ((1-rho)^2) ...
                + r2*ka*zeta*exp(ka/(1-rho)) / ((1-rho)^2) ...
                - r2*kc*zeta*exp(ka/(1-rho))/(1-rho)/rho ...
                - r2*kc*zeta/(rho^2)*igamma(0,-ka/(1-rho));
        end
        delta_k = - G_k; % Search direction
        %         save_rho = [save_rho, rho];
        k = k+1; % Increase iteration step
        rho = rho + step_k*delta_k; % Update rho
        disp(['k = ' num2str(k)]) % Current iteration step
    end
    rho_GRADIENT(ss) = rho;
    
end
disp(['k_end = ' num2str(k)]) % Final iteration step
disp(['rho_optimal_GRADIENT = ' num2str(rho)])
disp(['sum-throughput_GRADIENT = ' num2str(f_rho)])

% save rho.dat save_rho -ascii
