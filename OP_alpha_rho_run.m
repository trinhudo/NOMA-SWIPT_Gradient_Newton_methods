close all
clear all
%
%% Simulation parameters
%
K=3;
rho     = 0:0.2:1; % power splitting ratio
alpha   = 0:0.2:1; % time fraction for EH
PS_dB   = 30; % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = .5;
ncN     = .5;
naF     = .5;
ncF     = .5;
epsilon = 3; % pathloss exponent
%
dSF     = 1;
dSN     = .2;
dNF     = dSF - dSN;
%
lSF     = dSF^-epsilon;
lSN     = dSN^-epsilon;
lNF     = dNF^-epsilon;
%
eta     = 0.7; % energy conversion coefficient
pN      = 0.1; % power allocation coefficient
pF      = 1 - pN;
RthN    = 1; % bits/s/Hz
RthF    = ;
%
SimTimes = 10^0; % Monte-Carlo repetitions
%
%% Simulation
%
for xx = 1:length(alpha)
    disp(strcat('alpha=',num2str(alpha(xx))));
    for yy = 1:length(rho)
        disp(strcat('rho=',num2str(rho(yy))));
        %
        gthN = 2^(2*RthN/(1-alpha(xx))) - 1; % gamma_1
        gthF = 2^(RthF*2/(1-alpha(xx))) - 1; % gamma_2
        % channel modelling
        for ii = 1:K
            hSiF(:,ii) = sqrt(lSF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            hSiN(:,ii) = sqrt(lSN/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        end
        hNF = sqrt(lNF/2)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        % channel gains
        gSiN = abs(hSiN.^2);
        gSiF = abs(hSiF.^2);
        gNF  = abs(hNF.^2);
        % find the best far
        [gSsF,I] = max(gSiF,[],2);
        for tt = 1:SimTimes
            gSsN(tt,1) = gSiN(tt,I(tt));
        end
        % SNR modelling
        snrSsN_xF = (1-rho(yy)).*pF.*PS.*gSsN./...
            ((1-rho(yy)).*pN.*PS.*gSsN ...
            + (1-rho(yy))*naN + ncN);
        %
        snrSsN_xN = (1-rho(yy)).*pN.*PS.*gSsN/...
            (1-rho(yy))*naN + ncN;
        %
        snrSsF = pF.*PS.*gSsF./(pN.*PS.*gSsF + naF + ncF);
        %
        snrNF = eta.*PS.*gSsN.*gNF.*...
            (2*alpha(xx)/(1-alpha(xx))+rho(yy))/(naF + ncF);
        % count outage events
        count = 0;
        %
        for zz = 1:SimTimes
            %% for DF only
            if (snrSsN_xF(zz) >= gthF) && ...
                    (max(snrSsF(zz),snrNF(zz)) < gthF)
                count = count + 1;
            elseif (snrSsN_xF(zz) < gthF) && (snrSsF(zz) < gthF)
                count = count + 1;
            end
        end
        Pout_sim(xx,yy) = count/SimTimes;
        %             T_sim(aa,rr) = (1-Pout_sim(aa,rr))*RthF;
        %% Analytical Results
        a1 = (1-rho(yy))*pF*PS/((1-rho(yy))*naN + ncN);
        a2 = (1-rho(yy))*pN*PS/((1-rho(yy))*naN + ncN);
        b1 = pF * PS / (naF + ncF);
        b2 = pN * PS / (naF + ncF);
        c  = eta*PS*(2*alpha(xx)/(1-alpha(xx))+rho(yy))/(naF + ncF);
        %
        Phi1 = 1 - exp(-gthF/lSN/(a1-a2*gthF));
        %
        Phi2 = 0;
        %
        for jj = 0:K
            Phi2_temp = nchoosek(K,jj)*((-1)^jj)*...
                exp(-jj*gthF/lSF/(b1-b2*gthF));
            Phi2 = Phi2 + Phi2_temp;
        end
        %
        Theta1 = gthF/lSN/lNF/c*igamma(0,gthF/lSN/(a1-a2*gthF));
        % New approx. (2016-08-20)
        mu = gthF/(a1-a2*gthF);
        nu = 1/lSN;
        xi = gthF/lNF/c;
        %
        a = mu;
        b = xi;
        c = nu;
        %
        A1 = exp(-a*c)/c;
        A2 = -b*igamma(0,a*c);
        A3 = 0;
        for nn=2:2
            B1 = ((-1)^nn)*(b^nn)/(factorial(nn));
            B21 = exp(-(a*c));
            B22 = 0;
            for vv = 1:(nn-1)
                temp = (factorial(vv-1))*((-c)^(nn-vv-1))/...
                    ((factorial(nn-1))*(a^vv));
                B22 = B22 + temp;
            end
            B23 = ((-c)^(nn-1))/(factorial(nn-1))*(ei(-a*c));
            B2 = B21*B22-B23;
            A3 = A3 + B1*B2;
        end
        NewApprox = A3;
        %
        Pout_ana(xx,yy) = Phi2*(Phi1 + Theta1 - 1/lSN*NewApprox);
    end
end

for xx = 1:length(alpha)
    for yy = 1:length(rho)
        if (0 == isreal(Pout_ana(xx,yy)))
            Pout_ana(xx,yy) = 1;
        end
    end
end
%% plot
% load('dataSim')
surf(alpha,rho,Pout_ana,'facecolor','none','linestyle',':','marker','none')