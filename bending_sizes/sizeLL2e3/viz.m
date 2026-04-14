%% plot_dd_output_iter3.m
% Visualize regular history, nucleation-induced Ktip drops, peak envelope,
% cumulative free-surface slip-out count, and dislocation distribution evolution.
% Figure 3:
%   x-axis  = dislocation position
%   y-axis  = time
%   color   = dislocation ID
%   neutral axis marked by a vertical dashed line

clear; clc;

useAutoFind = true;
mainCsv = './output/outputKD_tau500.csv';
peakCsv = './output/KtipPeaks_tau500.csv';
disDir  = './disConfig';
inputFile = './input.txt';

saveFigures = false;
outDirFig = './figures';

% --------------------- general options ----------------------------------
plotDislocationFigure = true;
markerSizeDis = 28;
useAthermalShape = false;   % false: all circles; true: thermal=o, athermal=s

% --------------------- neutral-axis options -----------------------------
showNeutralAxis = true;

% If true, use the value below directly.
% If false, read input.txt and compute:
%   neutralAxisPos = crack_tip + bending_neutralXi * LL / abs(cos(theta))
useManualNeutralAxis = false;
neutralAxisPos = 0.0;

% default value if bending_neutralXi is absent in input.txt
default_bending_neutralXi = 0.5;
% ------------------------------------------------------------------------

if useAutoFind
    fileList = dir(fullfile('.', 'output', 'outputKD_tau*.csv'));
    if isempty(fileList)
        error('Main csv not found.');
    end
    [~, idx] = max([fileList.datenum]);
    mainCsv = fullfile(fileList(idx).folder, fileList(idx).name);

    peakList = dir(fullfile('.', 'output', 'KtipPeaks_tau*.csv'));
    if ~isempty(peakList)
        [~, idx2] = max([peakList.datenum]);
        peakCsv = fullfile(peakList(idx2).folder, peakList(idx2).name);
    end
end

% --------------------- read main table ----------------------------------
T = readtable(mainCsv);
T = sortrows(T, {'time','kInc','recordType'});

isRegular      = (T.recordType == 0) | (T.recordType == 3);  % for history plots
isRegularOut   = (T.recordType == 0);                        % only regular rows match evl_*.csv
isNucBefore    = (T.recordType == 1);
isNucAfter     = (T.recordType == 2);

tr             = T.time(isRegular);
Kapp_r         = T.Kapp_MPam05(isRegular);
Ktip_r         = T.Ktip_MPam05(isRegular);
Nd_r           = T.Nd(isRegular);
slipOutTotal_r = T.removedBoundaryTotal(isRegular);

Treg = T(isRegularOut, :);   % regular output rows only, for mapping disConfig snapshots

%% ==================== Figure 1 =========================================
fig1 = figure('Color','w','Position',[100 80 1000 720]);
hold on; box on;

plot(T.time, T.Ktip_MPam05, '-', 'LineWidth', 1.2, ...
    'DisplayName', 'K_{tip} full history');
plot(T.time, T.Kapp_MPam05, '--', 'LineWidth', 1.5, ...
    'DisplayName', 'K_{app}');
plot(T.time(isNucBefore), T.Ktip_MPam05(isNucBefore), 'o', ...
    'MarkerSize', 6, 'DisplayName', 'K_{tip} before nucleation');
plot(T.time(isNucAfter), T.Ktip_MPam05(isNucAfter), 's', ...
    'MarkerSize', 5, 'DisplayName', 'K_{tip} after nucleation');

if exist(peakCsv, 'file')
    P = readtable(peakCsv);
    if ~isempty(P)
        plot(P.time, P.Ktip_before_MPam05, '-', 'LineWidth', 2.0, ...
            'DisplayName', 'Peak envelope');
    end
end

xlabel('Time');
ylabel('SIF (MPa\cdotm^{0.5})');
title('K_{app}, K_{tip}, nucleation drops, and peak envelope');
legend('Location','best');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

%% ==================== Figure 2 =========================================
fig2 = figure('Color','w','Position',[120 90 1000 720]);

subplot(2,1,1); hold on; box on;
plot(tr, Nd_r, '-', 'LineWidth', 1.8);
xlabel('Time');
ylabel('N_d');
title('Dislocation number');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

subplot(2,1,2); hold on; box on;
plot(tr, slipOutTotal_r, '-', 'LineWidth', 1.8);
xlabel('Time');
ylabel('Cumulative slip-out count');
title('Dislocations slipped out of free surface');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

%% ==================== Figure 3 =========================================
% Dislocation distribution evolution:
% x-axis  = dislocation position
% y-axis  = time
% color   = dislocation ID
% neutral axis is marked by a vertical line

fig3 = [];

if plotDislocationFigure
    disFiles = dir(fullfile(disDir, 'evl_*.csv'));

    if isempty(disFiles)
        warning('No dislocation config files found in folder: %s. Figure 3 is skipped.', disDir);
    else
        % parse snapshot index from filenames
        snapIdx = nan(numel(disFiles), 1);
        for i = 1:numel(disFiles)
            snapIdx(i) = parseSnapshotIndex(disFiles(i).name);
        end

        validMask = ~isnan(snapIdx);
        disFiles  = disFiles(validMask);
        snapIdx   = snapIdx(validMask);

        [snapIdx, order] = sort(snapIdx);
        disFiles = disFiles(order);

        allTime = [];
        allPos  = [];
        allID   = [];
        allAth  = [];
        allVel  = [];
        allRss  = [];
        allSnap = [];

        for i = 1:numel(disFiles)
            outIdx = snapIdx(i);
            rowIdx = outIdx + 1;   % evl_0 <-> first regular row

            if rowIdx > height(Treg)
                warning('Snapshot %s has no matching regular row in main csv. Skipped.', disFiles(i).name);
                continue;
            end

            fpath = fullfile(disFiles(i).folder, disFiles(i).name);
            M = readmatrix(fpath);

            if isempty(M)
                continue;
            end

            % Ensure single-row csv is still treated as 1x5
            if isvector(M)
                M = reshape(M, 1, []);
            end

            if size(M,2) < 5
                warning('File %s has fewer than 5 columns. Skipped.', disFiles(i).name);
                continue;
            end

            nNow = size(M,1);

            allTime = [allTime; repmat(Treg.time(rowIdx), nNow, 1)];
            allID   = [allID;   M(:,1)];
            allAth  = [allAth;  M(:,2)];
            allPos  = [allPos;  M(:,3)];
            allVel  = [allVel;  M(:,4)];
            allRss  = [allRss;  M(:,5)];
            allSnap = [allSnap; repmat(outIdx, nNow, 1)];
        end

        if isempty(allTime)
            warning('No valid dislocation snapshot data found. Figure 3 is skipped.');
        else
            % Sort by time first, then by position
            A = [allTime, allPos, allID, allAth, allVel, allRss, allSnap];
            A = sortrows(A, [1 2]);

            allTime = A(:,1);
            allPos  = A(:,2);
            allID   = A(:,3);
            allAth  = A(:,4);
            allVel  = A(:,5);
            allRss  = A(:,6);
            allSnap = A(:,7);

            % --------- compute neutral-axis position ---------------------
            neutralAxisPos_used = NaN;
            if showNeutralAxis
                if useManualNeutralAxis
                    neutralAxisPos_used = neutralAxisPos;
                else
                    if exist(inputFile, 'file')
                        params = readInputTxtSimple(inputFile);

                        crack_tip = getFieldOrDefault(params, 'crack_tip', 0.0);
                        LL        = getFieldOrDefault(params, 'LL', NaN);
                        theta     = getFieldOrDefault(params, 'theta', pi/4);
                        bending_neutralXi = getFieldOrDefault(params, ...
                            'bending_neutralXi', default_bending_neutralXi);

                        if ~isnan(LL) && abs(cos(theta)) > 1e-12
                            slipLength_local = LL / abs(cos(theta));
                            neutralAxisPos_used = crack_tip + bending_neutralXi * slipLength_local;
                        end
                    end
                end
            end

            fig3 = figure('Color','w','Position',[150 110 1080 780]);
            hold on; box on;

            % Scatter: x = position, y = time
            if ~useAthermalShape
                scatter(allPos, allTime, markerSizeDis, allID, 'filled', ...
                    'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.92);
            else
                maskTherm  = (allAth == 0);
                maskAtherm = (allAth ~= 0);

                if any(maskTherm)
                    scatter(allPos(maskTherm), allTime(maskTherm), markerSizeDis, allID(maskTherm), ...
                        'o', 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.92);
                end
                if any(maskAtherm)
                    scatter(allPos(maskAtherm), allTime(maskAtherm), markerSizeDis+8, allID(maskAtherm), ...
                        's', 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.92);
                end
            end

            colormap(turbo(max(64, min(256, numel(unique(allID))))));
            cb = colorbar;
            cb.Label.String = 'Dislocation ID';

            % Mark neutral axis with a vertical dashed line
            if showNeutralAxis && ~isnan(neutralAxisPos_used)
                xline(neutralAxisPos_used, '--k', 'LineWidth', 1.8, ...
                    'Label', 'Neutral axis', ...
                    'LabelOrientation', 'horizontal', ...
                    'LabelVerticalAlignment', 'middle', ...
                    'HandleVisibility', 'off');
            end

            xlabel('Dislocation position');
            ylabel('Time');
            title('Dislocation distribution evolution (colored by ID)');
            set(gca, 'FontSize', 13, 'LineWidth', 1.1);

            xMin = min(allPos);
            xMax = max(allPos);
            yMin = min(allTime);
            yMax = max(allTime);

            dx = xMax - xMin;
            dy = yMax - yMin;
            if dx <= 0, dx = 1; end
            if dy <= 0, dy = 1; end

            xlim([xMin - 0.04*dx, xMax + 0.06*dx]);
            ylim([yMin - 0.02*dy, yMax + 0.03*dy]);

            if showNeutralAxis && ~isnan(neutralAxisPos_used)
                txt = sprintf(['Snapshots: %d   Total markers: %d   Max ID: %d\n' ...
                               'Neutral axis position = %.6g'], ...
                    numel(unique(allSnap)), numel(allID), max(allID), neutralAxisPos_used);
            else
                txt = sprintf('Snapshots: %d   Total markers: %d   Max ID: %d', ...
                    numel(unique(allSnap)), numel(allID), max(allID));
            end

            text(0.02, 0.98, txt, 'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 4);
        end
    end
end

%% ==================== save figures =====================================
if saveFigures
    if ~exist(outDirFig, 'dir')
        mkdir(outDirFig);
    end

    exportgraphics(fig1, fullfile(outDirFig, 'Ktip_drop_envelope.png'), 'Resolution', 300);
    exportgraphics(fig2, fullfile(outDirFig, 'Nd_slipout.png'), 'Resolution', 300);

    if ~isempty(fig3) && isgraphics(fig3)
        exportgraphics(fig3, fullfile(outDirFig, 'dislocation_distribution_vs_time_position.png'), 'Resolution', 300);
    end
end

%% ==================== local functions ==================================
function idx = parseSnapshotIndex(fname)
% parse evl_123.csv -> 123
    tok = regexp(fname, '^evl_(\d+)\.csv$', 'tokens', 'once');
    if isempty(tok)
        idx = NaN;
    else
        idx = str2double(tok{1});
    end
end

function params = readInputTxtSimple(fname)
% Read lines like:
%   key=value;
% and return a struct with numeric values when possible

    params = struct();

    fid = fopen(fname, 'r');
    if fid < 0
        warning('Cannot open input file: %s', fname);
        return;
    end

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end

        tline = strtrim(tline);
        if isempty(tline), continue; end

        % remove comments after %
        pct = strfind(tline, '%');
        if ~isempty(pct)
            tline = strtrim(tline(1:pct(1)-1));
        end
        if isempty(tline), continue; end

        tok = regexp(tline, '^\s*([A-Za-z_]\w*)\s*=\s*([^;]+)\s*;?\s*$', 'tokens', 'once');
        if isempty(tok), continue; end

        key = strtrim(tok{1});
        valStr = strtrim(tok{2});
        valNum = str2double(valStr);

        if ~isnan(valNum)
            params.(key) = valNum;
        else
            params.(key) = valStr;
        end
    end

    fclose(fid);
end

function val = getFieldOrDefault(S, fieldName, defaultVal)
    if isstruct(S) && isfield(S, fieldName)
        val = S.(fieldName);
    else
        val = defaultVal;
    end
end