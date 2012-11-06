function [t,f,median_f] = ppi_dyn(rawData,rawTime,bins,tFinal)
% Get dynamical data ready for plotting
% (c) G. Jack Peterson, 8/19/2011

optimax = 50;

csT = cumsum(rawTime')';
T = round(csT(:,50:50:length(csT))/bins)*bins;
% T = round(csT(:,1:length(csT))/bins)*bins;
maxT = max(max(T));
t = bins:bins:maxT;
tLen = length(t);
f = nan(optimax,tLen);

sizeT = size(T);

for ii = 1:sizeT(1)
    for jj = 1:tLen
%        % Median value of data associated with each t per simulation
% %        tCoords = T(ii,:) == t(jj);
% %        f(ii,jj) = median(nonzeros(rawData(ii,:).*tCoords(1:length(rawData(ii,:)))));
%        tempT = T(ii,:);
%        tempT(length(rawData(ii,:))+1:end) = [];
%        tCoords = find(tempT == t(jj));
%        f(ii,jj) = median(rawData(ii,tCoords));
        % Median value of data associated with each t per simulation
        tempT = T(ii,1:round(tFinal(ii)/50));
%         tempT(length(rawData(ii,:))+1:end) = [];
        tCoords = find(tempT == t(jj));
        if ~isempty(tCoords)
            f(ii,jj) = median(rawData(ii,tCoords));
        end

    end
end

% Compute median over simulations
median_f = nan(1,tLen);
for jj = 1:tLen
    temp = f(:,jj);
    if sum(~isnan(temp))
        temp(isnan(temp)) = [];
        median_f(jj) = median(temp);
    end
end
t(isnan(median_f)) = [];
median_f(isnan(median_f)) = [];