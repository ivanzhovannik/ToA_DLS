function ACF = TTTCU(dAT, dATminNs, dTauMinNs)

%% Time-Tag-To-Correlation-Unique Algorithm

% Given
% Vector of the arrival times AT = {t1,t2,...tN}
% For a given lag time Tau a vector is generated: AT_ = AT + Tau
% The second vector of the arrival times AT_ = {t1+Tau,t2+Tau,...tN+Tau}
% The algorithm makes the index of AT runs from the beginning until
% AT_(1)<=AT(j), if AT(j)=AT_(1) then ACF(Tau)++. Then the index of AT_
% runs from 1 to j1: AT_(j1)>=AT(j), if AT(j)=AT_(j1) then ACF++.

% k is the Tau index
% dTauMinNs is the minimum sampling time of the Tau vector in Nsec
% dATminNs is the sampling time of the Time vector in Nsec

ATLength = length(dAT);
AT = zeros([1 ATLength]); % defining the AT vector
AT(1) = dAT(1);
for n = 2:ATLength
    AT(n) = AT(n-1)+dAT(n);
end;
Tau = GenerateLagTimeVect(max(AT),2,dTauMinNs, dATminNs);
TauLength = length(Tau);
AT_ = zeros([1 ATLength]); % Default AT shifted vector
ACF_ = zeros([1 TauLength]); % Default Auto-correlation function 

%% Algorithm

for k = 1:TauLength
   AT_(:) = AT(:)+1*Tau(k);
   C = [AT,AT_];
   ACF_ (k) = 0;
   UniqueAT = unique(C);
   ACF_ (k) = ATLength - .5*length(UniqueAT);
end
ACF = [Tau'; ACF_']; % ./ATLength;

end