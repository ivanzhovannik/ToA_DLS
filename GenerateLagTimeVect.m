%% GenerateLagTimeVect algorithm

% Tmax is the length of the data vector, B ia an integer base number
% regulating length of the Tau vector generated with the algorithm
% dTauMinNs is the minimum sampling time of the Tau vector in Nsec
% dATminNs is the sampling time of the Time vector in Nsec

function [Tau] = GenerateLagTimeVect(Tmax,B,dTauMinNs,dATminNs)

%% Parameters entering

% Tmax = upper correlation time
TauRawLim = 500; % maximum Tau index
%B =  integer base number

%% Multi-tau definition

j = 1; % there was one
TauRaw = zeros([1 TauRawLim]);
N = dTauMinNs/dATminNs; % must be an integer value
TauRaw(1) = N;
TmaxDivN = Tmax/N;
while (TauRaw(j) < TmaxDivN) % Tau & Tmax have different sampling system!!!
j = j+1;
TauRaw(j) = TauRaw(j-1)+N*2^floor((j-1)/B);
end;
Tau (1:j) = TauRaw (1:j);

end