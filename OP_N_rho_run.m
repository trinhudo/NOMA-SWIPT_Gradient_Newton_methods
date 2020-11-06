close all
clear all
%% Simulation parameters
%
rho     = 0:.01:1;               % power splitting ratio
PS_dB   = 10;         % transmit SNR = Ps/N0 in dB
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
%
SimTimes = 10^0;            % Monte-Carlo repetitions
%
%% Simulation
%
for ss = 1:length(PS_dB)
    for rr = 1:length(rho)
        fprintf('rho = %d dB \n',rho(rr))
        % Channel modelling
        hSF = sqrt(lSF/2)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        hSN = sqrt(lSN/2)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        hNF = sqrt(lNF/2)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        % Channel gains
        gSN     = abs(hSN.^2);
        gSF     = abs(hSF.^2);
        gNF     = abs(hNF.^2);
        % SNRs
        snrSNxF = (1-rho(rr)).*pF.*PS(ss).*gSN./...
            ((1-rho(rr)).*pN.*PS(ss).*gSN ...
            + (1-rho(rr))*naN + ncN);
        snrSNxN = (1-rho(rr)).*pN.*PS(ss).*gSN/...
            ((1-rho(rr))*naN + ncN);
        snrSF   = pF.*PS(ss).*gSF./(pN.*PS(ss).*gSF + naF + ncF);
        snrNF   = eta.*PS(ss).*gSN.*gNF.*rho(rr)/(naF + ncF);
        %
        % Find the best antenna for User F based on end-to-end SNR
        snrFe2e = min(snrSNxF,max(snrNF,snrSF));
        % count outage events of User N
        count1 = snrSNxF < g2;
        count2 = snrSNxF >=g2 & snrSNxN < g1;
        OP_N_sim(rr) = sum(count1+count2)/SimTimes;
        %% Analysis
        a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN);
        a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
        b1 = pF * PS(ss) / (naF + ncF);
        b2 = pN * PS(ss) / (naF + ncF);
        c  = eta*rho(rr)*PS(ss)/(naF + ncF);
        mu_a = g2/(a1-a2*g2);
        mu_b = g2/(b1-b2*g2);
        theta = pF/pN;
        %
        if g2 < theta && mu_a >= g1/a2
            OP_N_ana(rr) = 1 - exp(-mu_a/lSN);
        elseif g2 < theta && mu_a < g1/a2
            OP_N_ana(rr) = 1 - exp(-g1/lSN/a2);
        elseif g2 >= theta
            OP_N_ana(rr) = 1;
        end
    end
end
%% plot
semilogy(rho,OP_N_sim,'o:',...
    rho,OP_N_ana,'*-')
xlabel('SNR (dB)')
ylabel('OP')

