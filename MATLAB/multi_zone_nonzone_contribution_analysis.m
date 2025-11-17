% ------------------------------------------------------------
% Multi-Zone Contribution Analysis (FAST + Interactive Component Selection)
% Author: Zahra Motahari Rad
% ------------------------------------------------------------

clear; clc;
filename = 'CAT-2Z-PP.txt';  % your input file

%% === 1. Read metadata ===
fid = fopen(filename,'r');
lines = textscan(fid,'%s','Delimiter','\n','Whitespace','');
lines = lines{1};
fclose(fid);

topEventName = strtrim(strsplit(lines{1}));
topEventName = topEventName{2};
topEventProb = strtrim(strsplit(lines{4}));
topEventProb = str2double(topEventProb{2});

fprintf('\nTop Event: %s\n', topEventName);
fprintf('Top Event Probability: %.8e\n\n', topEventProb);

%% === 2. Component list (fixed master list) ===
allComponents = { ...
    'B','MJ','FC','FS','PDB', ...                 % Zone 1 possible parts % 'B','MJ','FC','FS','PDB',
    'ESC1','ESC2','ESC3','ESC4','ESC5','ESC6', ... % ESCs 'ESC1','ESC2','ESC3','ESC4','ESC5','ESC6','ESC7','ESC8'
    'M1','M2','M3','M4','M5','M6', ...       % Motors 'M1','M2','M3','M4','M5','M6','M7','M8'
    'P1','P2','P3','P4','P5','P6'};          % Propellers 'P1','P2','P3','P4','P5','P6','P7','P8'

fprintf('Available components:\n');
for i = 1:numel(allComponents)
    fprintf('  %2d. %s\n', i, allComponents{i});
end
fprintf('\n');

%% === 3. Ask for zones and select components interactively ===
nZones = input('How many zones do you have? ');

zoneNames = cell(nZones,1);
zoneComponents = cell(nZones,1);

for i = 1:nZones
    zoneNames{i} = sprintf('Zone%d', i);
    fprintf('\nSelect components for %s:\n', zoneNames{i});
    fprintf('(Enter indices separated by spaces, e.g. [1 2 3 4])\n');
    idx = input('Indices: ');   % just one prompt now
    idx = unique(idx(idx>=1 & idx<=numel(allComponents)));
    comps = allComponents(idx);
    zoneComponents{i} = comps;
end

fprintf('\nZone configuration:\n');
for i = 1:nZones
    fprintf('  %s: %s\n', zoneNames{i}, strjoin(zoneComponents{i},', '));
end
fprintf('\n');


%% === 4. Read data lines (skip first 5 lines) ===
rawText = fileread(filename);
lines = regexp(rawText, '\r?\n', 'split');
lines = lines(6:end);
nLines = numel(lines);

%% === 5. Initialize totals ===
zoneProb = zeros(nZones,1);
zoneContrib = zeros(nZones,1);
zonebarProb = 0;
zonebarContrib = 0;

%% === 6. Process lines efficiently ===
tic;
for i = 1:nLines
    L = strtrim(lines{i});
    if isempty(L), continue; end
    parts = strsplit(L);
    if numel(parts) < 5, continue; end

    order = str2double(parts{2});
    prob = str2double(parts{3});
    contrib = str2double(parts{4});
    if isnan(prob) || isnan(contrib) || order <= 0, continue; end

    events = parts(5:end);
    nE = numel(events);
    shareP = prob / nE;
    shareC = contrib / nE;

    for j = 1:nE
        evt = strtrim(events{j});
        dotIdx = strfind(evt,'.');
        if isempty(dotIdx)
            comp = evt; mode = '';
        else
            comp = evt(1:dotIdx(1)-1);
            mode = evt(dotIdx(1)+1:end);
        end

        isZoneEvt = ~isempty(regexpi(mode,'(Temp|Vibra)$'));

        if isZoneEvt
            % Assign to the correct zone based on component
            assigned = false;
            for z = 1:nZones
                if any(strcmpi(comp, zoneComponents{z}))
                    zoneProb(z) = zoneProb(z) + shareP;
                    zoneContrib(z) = zoneContrib(z) + shareC;
                    assigned = true;
                    break;
                end
            end
            if ~assigned
                zonebarProb = zonebarProb + shareP;
                zonebarContrib = zonebarContrib + shareC;
            end
        else
            % Non-zone event
            zonebarProb = zonebarProb + shareP;
            zonebarContrib = zonebarContrib + shareC;
        end
    end
end
tElapsed = toc;
fprintf('Processing time: %.2f seconds for %d cutsets.\n', tElapsed, nLines);

%% === 7. Build summary ===
totalContrib = sum(zoneContrib) + zonebarContrib;

Category = [zoneNames; {'ZoneBar'}];
Total_Probability = [zoneProb; zonebarProb];
Total_Contribution = [zoneContrib; zonebarContrib];
Contribution_Pct = 100 * (Total_Contribution / totalContrib);
TopEventName = repmat({topEventName}, numel(Category),1);
TopEventProbability = repmat(topEventProb, numel(Category),1);

summaryTable = table(TopEventName, TopEventProbability, Category, ...
    Total_Probability, Total_Contribution, Contribution_Pct);

disp('================ SUMMARY ================');
disp(summaryTable);

%% === 8. Save results ===
safeEventName = regexprep(topEventName,'\W','_');
summaryFile = sprintf('ZoneContributionSummary-%s.csv',safeEventName);
writetable(summaryTable, summaryFile);

fprintf('\n✅ Multi-zone summary saved: %s\n', summaryFile);



%% === NEW: Component-Level Zone Contribution Ranking ===
componentMaps = cell(nZones,1);
for z = 1:nZones
    componentMaps{z} = containers.Map('KeyType','char','ValueType','double');
end

rawText = fileread(filename);
linesEvt = regexp(rawText, '\r?\n', 'split');
linesEvt = linesEvt(6:end);
nLinesEvt = numel(linesEvt);

for i = 1:nLinesEvt
    L = strtrim(linesEvt{i});
    if isempty(L), continue; end
    parts = strsplit(L);
    if numel(parts) < 5, continue; end

    contrib = str2double(parts{4});
    if isnan(contrib), continue; end

    events = parts(5:end);
    nE = numel(events);
    shareC = contrib / nE;

    for k = 1:nE
        ev = events{k};

        % Zone event check
        if isempty(regexp(ev,'(Temp|Vibra)$','once'))
            continue;
        end

        % Extract component name
        comp = regexp(ev, '^[^.]+', 'match', 'once');

        % Find which zone this component belongs to
        for z = 1:nZones
            if any(strcmpi(comp, zoneComponents{z}))
                if componentMaps{z}.isKey(comp)
                    componentMaps{z}(comp) = componentMaps{z}(comp) + shareC;
                else
                    componentMaps{z}(comp) = shareC;
                end
                break;
            end
        end
    end
end

%% === Export each zone table ===
safeEventName = regexprep(topEventName,'\W','_');
outFile = sprintf('ZoneContribution-%s.xlsx', safeEventName);

for z = 1:nZones
    cmap = componentMaps{z};
    comps = cmap.keys;
    vals  = cmap.values;

    if isempty(comps), continue; end

    compTable = table(comps', cell2mat(vals'), ...
        'VariableNames', {'Component','Contribution'});

    % Sort
    compTable = sortrows(compTable, 'Contribution', 'descend');

    % Percent contributions
    compTable.PctOfZone = 100 * compTable.Contribution / zoneContrib(z);
    compTable.PctOfAll  = 100 * compTable.Contribution / totalContrib;

    % Save
    writetable(compTable, outFile, 'Sheet', sprintf('%s_Components', zoneNames{z}));

    fprintf('✅ Component ranking saved for %s\n', zoneNames{z});
end

fprintf('\n✅ Component-level tables added to: %s\n\n', outFile);



%% === NEW: Component-Level NON-Zone Contribution Ranking ===
zoneBarMap = containers.Map('KeyType','char','ValueType','double');

for i = 1:nLinesEvt
    L = strtrim(linesEvt{i});
    if isempty(L), continue; end
    parts = strsplit(L);
    if numel(parts) < 5, continue; end

    contrib = str2double(parts{4});
    if isnan(contrib), continue; end

    events = parts(5:end);
    nE = numel(events);
    shareC = contrib / nE;

    for k = 1:nE
        ev = strtrim(events{k});

        % if event IS NOT zone-related
        if isempty(regexp(ev,'(Temp|Vibra)$','once'))

            % component name
            comp = regexp(ev, '^[^.]+', 'match', 'once');

            % register contribution
            if zoneBarMap.isKey(comp)
                zoneBarMap(comp) = zoneBarMap(comp) + shareC;
            else
                zoneBarMap(comp) = shareC;
            end
        end
    end
end

% === Export Non-Zone Table ===
if zoneBarMap.Count > 0
    comps = zoneBarMap.keys;
    vals  = zoneBarMap.values;

    nonZoneTable = table(comps', cell2mat(vals'), ...
        'VariableNames', {'Component','Contribution'});

    % Sort
    nonZoneTable = sortrows(nonZoneTable, 'Contribution', 'descend');

    % Percent contributions
    nonZoneTable.PctOfNonZone = 100 * nonZoneTable.Contribution / zonebarContrib;
    nonZoneTable.PctOfAll     = 100 * nonZoneTable.Contribution / totalContrib;

    % Save
    writetable(nonZoneTable, outFile, 'Sheet', 'ZoneBar_Components');

    fprintf('✅ Component ranking saved for NON-zone cutsets (ZoneBar)\n');
else
    fprintf('⚠ No non-zone cutsets found.\n');
end

