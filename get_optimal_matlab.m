close all
clear all
%% Simulation parameters
%
rho     = 0:.01:1;          % power splitting ratio
PS_dB   = 10;               % transmit SNR = Ps/N0 in dB
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
RthF    = 1;               % target data rate of User N bits/s/Hz
g1 = 2^(2*RthN) - 1; % gamma_2 for User F
g2 = 2^(2*RthF) - 1; % gamma_2 for User F
%
%% Simulation
%
for ss = 1:length(PS_dB)
    for rr = 1:length(rho)
        fprintf('rho = %1f\n',rho(rr))
        %% Analysis
        a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN); % omitting (1-rho) at noise power
        a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
        b1 = pF * PS(ss) / (naF + ncF);
        b2 = pN * PS(ss) / (naF + ncF);
        c  = eta*rho(rr)*PS(ss)/(naF + ncF);
        mu_a = g2/(a1-a2*g2);
        mu_b = g2/(b1-b2*g2);
        theta = pF/pN;
        %
        OP_F = 1 - exp(-mu_a/lSN - mu_b/lSF)...
            - (1 - exp(-mu_b/lSF)) ...
            * (exp(-mu_a/lSN) - g2/lSN/lNF/c*igamma(0,mu_a/lSN));
        throughput_F(rr) = (1-OP_F)*RthF;
        %
        if g2/(pF-pN*g2) >= g1/pN % omitting the condition of g2<theta
            OP_N = 1 - exp(-mu_a/lSN);
            T_sum_ana(rr) = (1-OP_N)*RthN + (1-OP_F)*RthF;
        elseif g2/(pF-pN*g2) < g1/pN
            OP_N = 1 - exp(-g1/lSN/a2);
            T_sum_ana(rr) = (1-OP_F)*RthF + (1-OP_N)*RthN;
        elseif g2 >= theta
            T_sum_ana(rr) = 2;
        end
    end
    ID_rho = find(T_sum_ana == max(T_sum_ana));
    rho_opt = rho(ID_rho);
    sum_throughput_optimal(ss) = T_sum_ana(ID_rho);
    throughput_F_opt(ss) = throughput_F(ID_rho);
end


%% plot
plot(rho,T_sum_ana,'b-')
hold on
plot(rho_opt,T_sum_ana(ID_rho),'p')
%
% plot(PS_dB,sum_throughput_optimal,'s-')
% xlabel('P_S (dB)')
% ylabel('Optimal sum-throughput')
% axis([0 1 0 1])

