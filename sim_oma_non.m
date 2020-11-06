close all
clear all
%% Simulation parameters
%
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
%
g1_oma = 2^(2*RthN) - 1;
g1_non = 2^RthN - 1;
g2_oma  = 2^(2*RthF) - 1;
g2_non  = 2^RthF - 1;
SimTimes = 1e6;
%% Simulation
%
for ss = 1:length(PS_dB)
    fprintf('SNR = %d dB \n',PS_dB(ss))
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
    % OMA
    snrSN_oma = PS(ss).*gSN./(naN + ncN);
    snrSF_oma = PS(ss).*gSF./(naF + ncF);
    % non-coop
    snrSNxF = pF.*PS(ss).*gSN./(pN.*PS(ss).*gSN + naN + ncN);
    snrSNxN = pN.*PS(ss).*gSN/ (naN + ncN);
    snrSF   = pF.*PS(ss).*gSF./(pN.*PS(ss).*gSF + naF + ncF);
    % count outage events
    count_N_oma = snrSN_oma < g1_oma;
    OP_N_oma(ss) = sum(count_N_oma)/SimTimes;
    %
    count_F_oma = snrSF_oma < g2_oma;
    OP_F_oma(ss) = sum(count_F_oma)/SimTimes;
    %
    count_N_non1 = snrSNxF < g2_non;
    count_N_non2 = snrSNxF >= g2_non & snrSNxN < g1_non;
    OP_N_non(ss) = sum(count_N_non1+count_N_non2)/SimTimes;
    %
    count_F_non = min(snrSNxF,snrSF) < g2_non;
    OP_F_non(ss) = sum(count_F_non)/SimTimes;
%     % sum througput
%     sum_oma(ss) = (1-OP_N_oma)*RthN + (1-OP_F_oma)*RthF;
%     sum_non(ss) = (1-OP_N_non)*RthN + (1-OP_F_non)*RthF;
    %
    N_oma(ss) = (1-OP_N_oma(ss))*RthN;
    N_non(ss) = (1-OP_N_non(ss))*RthN;
    %
    F_oma(ss) = (1-OP_F_oma(ss))*RthF;
    F_non(ss) = (1-OP_F_non(ss))*RthF;
end
% %% plot
% plot(PS_dB,sum_oma,'o-',...
%     PS_dB,sum_non,'*-')
% hold on 
% plot(PS_dB,F_oma,'o-',...
%     PS_dB,F_non,'*-')
% legend('OMA','Non-coop')
% xlabel('P_S (dB)')
% ylabel('sum-throughput')