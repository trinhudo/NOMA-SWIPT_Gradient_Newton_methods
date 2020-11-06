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
        %
        if g2/(pF-pN*g2) >= g1/pN % omitting the condition of g2<theta
            OP_N = 1 - exp(-mu_a/lSN);
            T_sum_ana(rr) = (1-OP_N)*RthN + (1-OP_F)*RthF;
            obj_fun(rr) = r1*exp(nu1/(1-rho(rr))) ...
                + r2*exp(ka/(1-rho(rr))+kb) ...
                + r2*zeta*exp(ka/(1-rho(rr))) ...
                + r2*zeta*kc/rho(rr)*igamma(0,-ka/(1-rho(rr)));
        elseif g2/(pF-pN*g2) < g1/pN
            OP_N = 1 - exp(-g1/lSN/a2);
            T_sum_ana(rr) = (1-OP_F)*RthF + (1-OP_N)*RthN;
            obj_fun(rr) = r1*exp(nu2/(1-rho(rr))) ...
                + r2*exp(ka/(1-rho(rr))+kb) ...
                + r2*zeta*exp(ka/(1-rho(rr))) ...
                + r2*zeta*kc/rho(rr)*igamma(0,-ka/(1-rho(rr)));
        elseif g2 >= theta
            T_sum_ana(rr) = 2;
        end
    end
end
% ID_rho = find(T_sum_ana == max(T_sum_ana));
%         rho_opt = rho(ID_rho);
rho_gradient = 0.2117;
T_gradient = 1.8447;
%% plot
plot(rho,T_sum_ana,'b-', rho, -obj_fun,'r--')
hold on
plot(rho_gradient,T_gradient,'p')
%
legend('Sum-throughput','Ojective function','Optimal value', ...
    'location','southwest')
xlabel('\rho')
ylabel('Sum-throughput (bits/s/Hz)')
axis([0 1 0 2])
xticks(0:0.2:1)
yticks(0:.5:2)

