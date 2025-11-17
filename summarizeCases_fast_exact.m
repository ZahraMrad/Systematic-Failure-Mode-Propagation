%% summarizeCases_fast_exact.m
% -------------------------------------------------------------------------
% Exact output match to the “pre” version, but MUCH faster.
% Multi-rotor controllability case summarization for SFMP / Chapter 3.
%
% This script:
%   - Loads large rotor-state scenario tables from Excel
%   - Compresses cases using fixed-point merging of bitmask states
%   - Exports minimized switch-case logic for AltaRica
%
% Author: Zahra Motahari Rad
% -------------------------------------------------------------------------

clearvars -except summarizeCases_fast_exact; clc;

%% === Config ===
filename   = 'Control-MATLAB-Coax.xlsx';   % your file
rotorCols  = 23:30;                        % 23:30 for 8 rotors (Coax), 25:30 for 6 rotor (Hexa)
outCol     = 31;                           % AE
TOK_ORDER  = ["CATCtrl","HZDCtrl","MJRCtrl","MNRCtrl","Working"];

% Token -> bit mask mapping
bitvals = uint8([1 2 4 8 16]);             % fixed order
tok2bit = containers.Map(TOK_ORDER, num2cell(bitvals));

%% === Load Excel ===
fprintf('Reading %s ...\n', filename);
T = readtable(filename, 'ReadVariableNames', false);
T = T(2:end, :);   % row 11 onward (your original logic)

Rall = height(T);
fprintf('Rows after trim: %d\n', Rall);

rotorTbl   = T(:, rotorCols);
outputsCol = string(table2array(T(:, outCol)));
clear T;

%% === Normalize rotor tokens (string -> uint8 bitmask) ===
C = numel(rotorCols);
R = height(rotorTbl);
masks = zeros(R, C, 'uint8');

for j = 1:C
    col = string(table2array(rotorTbl(:, j)));

    for k = 1:numel(TOK_ORDER)
        m = strcmpi(col, TOK_ORDER(k));
        masks(m, j) = bitor(masks(m, j), bitvals(k));
    end

    mUnknown = (masks(:, j) == 0);
    masks(mUnknown, j) = uint8(16); % Working
end

labels = strtrim(outputsCol);
labels(ismissing(labels) | labels=="") = "Working";
cats = unique(labels, 'stable');

fprintf('Categories found: %d\n', numel(cats));

%% === Setup parallel (optional)
try
    p = gcp('nocreate');
    if isempty(p)
        parpool('local');   % Worker processes (supports fopen)
    end
catch
    warning('Could not start parallel pool. Running serial.');
end

%% === Process each category (can be parallel)
for kk = 1:numel(cats)
    label = cats(kk);
    idx = strcmp(labels, label);
    M = masks(idx, :);

    if isempty(M)
        continue;
    end

    fprintf('[%s] Start (%d rows)\n', label, size(M,1));

    %% === FIXED-POINT MERGING ===
    changed = true;
    while changed
        changed = false;

        for c = 1:C
            % Build keys for all-but-column-c
            if c == 1
                other = M(:, 2:end);
            elseif c == C
                other = M(:, 1:C-1);
            else
                other = [M(:,1:c-1) M(:,c+1:C)];
            end

            key = join(string(other), ',');
            key = join(key, ',', 2);

            [G, grpIdx] = findgroups(key);
            mc = M(:, c);

            counts = histcounts(G, 0.5:1:(max(G)+0.5));
            multi = find(counts > 1);

            if isempty(multi)
                continue;
            end

            changed = true;
            keep = true(size(M,1),1);
            newRows = zeros(numel(multi), C, 'uint8');
            w = 0;

            for g = multi
                rix = find(G == g);

                mergedC = uint8(0);
                for t = 1:numel(rix)
                    mergedC = bitor(mergedC, mc(rix(t)));
                end

                rep = M(rix(1), :);
                rep(c) = mergedC;

                w = w + 1;
                newRows(w,:) = rep;
                keep(rix) = false;
            end

            M = [M(keep,:); newRows(1:w,:)];

            % Deduplicate
            keyFull = join(string(M), ',');
            keyFull = join(keyFull, ',', 2);
            [~, ia] = unique(keyFull, 'stable');
            M = M(ia,:);
        end
    end

    %% === Emit AltaRica code ===
    nRows = size(M,1);
    lines = strings(nRows,1);

    for r = 1:nRows
        terms = strings(1,C);
        for j = 1:C
            mj = M(r,j);
            tokset = bitmaskToTokens(mj, TOK_ORDER, bitvals);
            if numel(tokset) == 1
                terms(j) = sprintf('(C_inFlowM%d == %s)', j, tokset);
            else
                parts = compose('(C_inFlowM%d == %s)', j, tokset);
                terms(j) = "(" + strjoin(parts, " or ") + ")";
            end
        end
        lines(r) = "(" + strjoin(terms, " and ") + ")";
    end

    allCases = strjoin(lines, " or ");
    code = sprintf(['C_outFlowC := switch {\n' ...
                    '    case %s : %s\n' ...
                    '    default : Working\n};'], ...
                    allCases, label);

    fname = sprintf('Switch_%s_summarized.txt', label);
    fid = fopen(fname,'w');
    fprintf(fid, '%s\n', code);
    fclose(fid);

    fprintf('[%s] Done. Merged %d rows.\n', label, nRows);
end

disp('✅ ALL DONE ✅');

%% ------------------------------------------------------------------------
%% Helper: decode bitmask to tokens
%% ------------------------------------------------------------------------
function toks = bitmaskToTokens(m, TOK_ORDER, bitvals)
    mask = (bitand(m, bitvals) ~= 0);
    toks = TOK_ORDER(mask);
end
