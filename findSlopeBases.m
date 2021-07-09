function slopeBases = findSlopeBases(x,peaksIn)

% Find positive going signal
y = [diff(x); 0];
y(y<0) = 0;
y(y>0) = 1;
y([1 numel(y)]) = 0;
candidates = strfind(y', [0 1])';

slopeBases = [];
for i = 1:numel(peaksIn)
    currentCandidates = candidates((peaksIn(i)-candidates)>0);
    [~,idx] = min(peaksIn(i)-currentCandidates);
    
    if ~isempty(idx) & sum(slopeBases==currentCandidates(idx))==0
        slopeBases(i) = peaksIn(i) - ...
            (peaksIn(i)-currentCandidates(idx));
    else % no slopeBase found, or is already taken
        slopeBases(i) = NaN; % infer this value later by taking median of successful pairings
    end
end

% Infer any missing slopeBases
idx = isnan(slopeBases);
if ~isempty(idx)
    slopeBases(idx) = (peaksIn(idx)-nanmedian(peaksIn-slopeBases'));
end

end