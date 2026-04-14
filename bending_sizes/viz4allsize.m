%% plot_batch_fracture_vs_LL_with_envelope.m
% 放在主文件夹下运行
% 自动读取各 sizeLLxe3 子文件夹中的 output/outputKD_tau*.csv
% 提取断裂时刻的 Kapp 和位错数量 N(Nd)
% 绘制线性图、双log图
% 并新增：每个算例的 Kapp-Ktip 上包络线子图对比

clear; clc; close all;

%% ================= 用户可调参数 =================
rootDir = pwd;                 % 主文件夹
saveFigures = true;            % 是否保存图片
saveSummaryTable = true;       % 是否保存汇总表
figDir = fullfile(rootDir, 'figures_batch');

% 若自动搜索不到，可手动指定 case 名称
useManualCaseList = false;
manualCaseList = { ...
    'sizeLL2e3', ...
    'sizeLL8e3', ...
    'sizeLL32e3', ...
    'sizeLL64e3'};

% 双log图是否加幂律拟合
showFitOnLoglog = true;

% 拟合时是否只用正值点
fitRequirePositive = true;

% ===== 新增：Kapp-Ktip 包络线参数 =====
makeEnvelopeFigure = true;     % 是否绘制 Kapp-Ktip 上包络线对比图
saveEnvelopeData   = true;     % 是否保存每个算例的包络线数据

% 包络线提取方式：
% 'binmax' : 按 Kapp 分箱，每箱取最大 Ktip（推荐）
envelopeMethod = 'binmax';

nEnvelopeBins = 80;            % 分箱数量（binmax 方法）
minPointsPerBin = 1;           % 每箱最少点数
smoothEnvelope = false;        % 是否对包络线做轻微平滑
smoothWindow    = 5;           % 平滑窗口（点数）

% 原始散点是否显示
showRawScatterInEnvelope = true;

% 子图布局
nColsEnvelope = 2;
% ==============================================

if saveFigures && ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% ================= 搜索算例文件夹 =================
if useManualCaseList
    caseDirs = {};
    for i = 1:numel(manualCaseList)
        cdir = fullfile(rootDir, manualCaseList{i});
        if exist(cdir, 'dir')
            caseDirs{end+1} = cdir; %#ok<SAGROW>
        else
            warning('Folder not found: %s', cdir);
        end
    end
else
    D = dir(rootDir);
    caseDirs = {};
    for i = 1:numel(D)
        if ~D(i).isdir
            continue;
        end
        namei = D(i).name;
        if strcmp(namei,'.') || strcmp(namei,'..')
            continue;
        end

        if ~isempty(regexp(namei, '^sizeLL\d+e3$', 'once')) || ...
           ~isempty(regexp(namei, '^size\d+e3$', 'once'))
            caseDirs{end+1} = fullfile(rootDir, namei); %#ok<SAGROW>
        end
    end
end

if isempty(caseDirs)
    error('未找到任何 sizeLLxe3 或 sizexe3 子文件夹。');
end

%% ================= 读取每个算例 =================
result = struct( ...
    'caseName', {}, ...
    'LL', {}, ...
    'csvPath', {}, ...
    'fractureTime', {}, ...
    'fractureKinc', {}, ...
    'recordType', {}, ...
    'Kapp_MPam05', {}, ...
    'Nd', {}, ...
    'KappHist', {}, ...
    'KtipHist', {}, ...
    'envX', {}, ...
    'envY', {} );

for i = 1:numel(caseDirs)
    caseDir = caseDirs{i};
    [~, caseName] = fileparts(caseDir);

    % 从文件夹名解析 LL
    LL = parseLLfromFolderName(caseName);
    if isnan(LL)
        warning('无法从文件夹名解析 LL，跳过: %s', caseName);
        continue;
    end

    % 优先找 output/outputKD_tau*.csv
    csvList = dir(fullfile(caseDir, 'output', 'outputKD_tau*.csv'));

    % 兼容：如果 output 子目录里没有，再试根目录
    if isempty(csvList)
        csvList = dir(fullfile(caseDir, 'outputKD_tau*.csv'));
    end

    if isempty(csvList)
        warning('未找到 outputKD_tau*.csv : %s', caseName);
        continue;
    end

    % 取最新文件
    [~, idxLatest] = max([csvList.datenum]);
    csvPath = fullfile(csvList(idxLatest).folder, csvList(idxLatest).name);

    T = readtable(csvPath);
    if isempty(T)
        warning('CSV 为空: %s', csvPath);
        continue;
    end

    % 必要列检查
    requiredCols = {'Kapp_MPam05','Nd'};
    for k = 1:numel(requiredCols)
        if ~ismember(requiredCols{k}, T.Properties.VariableNames)
            error('文件 %s 缺少列: %s', csvPath, requiredCols{k});
        end
    end

    % 寻找 Ktip 列（做兼容）
    ktipColName = findKtipColumn(T.Properties.VariableNames);
    if isempty(ktipColName)
        warning('文件 %s 中未找到 Ktip 列，包络线图将跳过该算例。', csvPath);
    end

    % 断裂时刻：优先 final row(recordType==3)，否则取最后一行
    if ismember('recordType', T.Properties.VariableNames)
        idxFinal = find(T.recordType == 3);
    else
        idxFinal = [];
    end

    if ~isempty(idxFinal)
        idxFrac = idxFinal(end);
    else
        idxFrac = height(T);
    end

    % 若最后行存在 NaN，则向前找最近的有效行
    while idxFrac >= 1
        if ~isnan(T.Kapp_MPam05(idxFrac)) && ~isnan(T.Nd(idxFrac))
            break;
        end
        idxFrac = idxFrac - 1;
    end

    if idxFrac < 1
        warning('文件中没有有效断裂行: %s', csvPath);
        continue;
    end

    % 历史 Kapp-Ktip
    KappHist = [];
    KtipHist = [];
    envX = [];
    envY = [];

    if ~isempty(ktipColName)
        KappHist = T.Kapp_MPam05;
        KtipHist = T.(ktipColName);

        validHist = isfinite(KappHist) & isfinite(KtipHist);
        KappHist = KappHist(validHist);
        KtipHist = KtipHist(validHist);

        if ~isempty(KappHist)
            [envX, envY] = computeUpperEnvelope(KappHist, KtipHist, ...
                envelopeMethod, nEnvelopeBins, minPointsPerBin, ...
                smoothEnvelope, smoothWindow);
        end
    end

    one.caseName = caseName;
    one.LL = LL;
    one.csvPath = csvPath;

    if ismember('time', T.Properties.VariableNames)
        one.fractureTime = T.time(idxFrac);
    else
        one.fractureTime = NaN;
    end

    if ismember('kInc', T.Properties.VariableNames)
        one.fractureKinc = T.kInc(idxFrac);
    else
        one.fractureKinc = NaN;
    end

    if ismember('recordType', T.Properties.VariableNames)
        one.recordType = T.recordType(idxFrac);
    else
        one.recordType = NaN;
    end

    one.Kapp_MPam05 = T.Kapp_MPam05(idxFrac);
    one.Nd = T.Nd(idxFrac);

    one.KappHist = KappHist;
    one.KtipHist = KtipHist;
    one.envX = envX;
    one.envY = envY;

    result(end+1) = one; %#ok<SAGROW>
end

if isempty(result)
    error('没有成功读取到任何有效算例。');
end

%% ================= 整理结果表 =================
caseName_col = string({result.caseName})';
LL_col       = [result.LL]';
time_col     = [result.fractureTime]';
kInc_col     = [result.fractureKinc]';
rtype_col    = [result.recordType]';
Kapp_col     = [result.Kapp_MPam05]';
Nd_col       = [result.Nd]';
csv_col      = string({result.csvPath})';

Summary = table(caseName_col, LL_col, time_col, kInc_col, rtype_col, Kapp_col, Nd_col, csv_col, ...
    'VariableNames', {'caseName','LL','fractureTime','fractureKinc','recordType','Kapp_MPam05','Nd','csvPath'});

Summary = sortrows(Summary, 'LL');

disp('===== Fracture summary =====');
disp(Summary);

if saveSummaryTable
    outCsv = fullfile(rootDir, 'summary_fracture_Kapp_N_vs_LL.csv');
    writetable(Summary, outCsv);
    fprintf('汇总表已保存: %s\n', outCsv);
end

%% ================= Figure 1: 线性坐标 =================
fig1 = figure('Color','w','Position',[100 80 900 760]);

subplot(2,1,1); hold on; box on;
plot(Summary.LL, Summary.Kapp_MPam05, '-o', ...
    'LineWidth', 1.8, 'MarkerSize', 7);
xlabel('LL');
ylabel('K_{app} at fracture (MPa\cdotm^{0.5})');
title('Fracture K_{app} vs LL (linear scale)');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

for i = 1:height(Summary)
    text(Summary.LL(i), Summary.Kapp_MPam05(i), ['  ' char(Summary.caseName(i))], ...
        'FontSize', 10, 'Interpreter', 'none');
end

subplot(2,1,2); hold on; box on;
plot(Summary.LL, Summary.Nd, '-s', ...
    'LineWidth', 1.8, 'MarkerSize', 7);
xlabel('LL');
ylabel('N_d at fracture');
title('Fracture dislocation number N vs LL (linear scale)');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

for i = 1:height(Summary)
    text(Summary.LL(i), Summary.Nd(i), ['  ' char(Summary.caseName(i))], ...
        'FontSize', 10, 'Interpreter', 'none');
end

%% ================= Figure 2: 双log坐标 =================
fig2 = figure('Color','w','Position',[120 90 900 760]);

subplot(2,1,1); hold on; box on;
maskK = isfinite(Summary.LL) & isfinite(Summary.Kapp_MPam05);
if fitRequirePositive
    maskK = maskK & (Summary.LL > 0) & (Summary.Kapp_MPam05 > 0);
end

loglog(Summary.LL(maskK), Summary.Kapp_MPam05(maskK), 'o-', ...
    'LineWidth', 1.8, 'MarkerSize', 7);

xlabel('LL');
ylabel('K_{app} at fracture (MPa\cdotm^{0.5})');
title('Fracture K_{app} vs LL (log-log)');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

if showFitOnLoglog && nnz(maskK) >= 2
    x = Summary.LL(maskK);
    y = Summary.Kapp_MPam05(maskK);

    p = polyfit(log10(x), log10(y), 1);
    slopeK = p(1);
    interceptK = p(2);

    xfit = logspace(log10(min(x)), log10(max(x)), 200);
    yfit = 10^(interceptK) * xfit.^slopeK;
    loglog(xfit, yfit, '--', 'LineWidth', 1.5);

    txtK = sprintf('fit: K_{app} \\propto LL^{%.3f}', slopeK);
    text(0.05, 0.92, txtK, 'Units', 'normalized', 'FontSize', 11, ...
        'BackgroundColor', 'w', 'Margin', 3);
end

for i = 1:height(Summary)
    if maskK(i)
        text(Summary.LL(i), Summary.Kapp_MPam05(i), ['  ' char(Summary.caseName(i))], ...
            'FontSize', 10, 'Interpreter', 'none');
    end
end

subplot(2,1,2); hold on; box on;
maskN = isfinite(Summary.LL) & isfinite(Summary.Nd);
if fitRequirePositive
    maskN = maskN & (Summary.LL > 0) & (Summary.Nd > 0);
end

loglog(Summary.LL(maskN), Summary.Nd(maskN), 's-', ...
    'LineWidth', 1.8, 'MarkerSize', 7);

xlabel('LL');
ylabel('N_d at fracture');
title('Fracture dislocation number N vs LL (log-log)');
set(gca, 'FontSize', 13, 'LineWidth', 1.1);

if showFitOnLoglog && nnz(maskN) >= 2
    x = Summary.LL(maskN);
    y = Summary.Nd(maskN);

    p = polyfit(log10(x), log10(y), 1);
    slopeN = p(1);
    interceptN = p(2);

    xfit = logspace(log10(min(x)), log10(max(x)), 200);
    yfit = 10^(interceptN) * xfit.^slopeN;
    loglog(xfit, yfit, '--', 'LineWidth', 1.5);

    txtN = sprintf('fit: N \\propto LL^{%.3f}', slopeN);
    text(0.05, 0.92, txtN, 'Units', 'normalized', 'FontSize', 11, ...
        'BackgroundColor', 'w', 'Margin', 3);
end

for i = 1:height(Summary)
    if maskN(i)
        text(Summary.LL(i), Summary.Nd(i), ['  ' char(Summary.caseName(i))], ...
            'FontSize', 10, 'Interpreter', 'none');
    end
end

%% ================= Figure 3: Kapp-Ktip 上包络线子图 =================
if makeEnvelopeFigure
    validEnv = arrayfun(@(s) ~isempty(s.envX) && ~isempty(s.envY), result);
    resultEnv = result(validEnv);

    if isempty(resultEnv)
        warning('没有可用于绘制 Kapp-Ktip 包络线的算例。');
    else
        % 按 LL 排序，便于对子图顺序一致
        [~, orderEnv] = sort([resultEnv.LL]);
        resultEnv = resultEnv(orderEnv);

        % 统一坐标范围
        allKapp = [];
        allKtip = [];
        for i = 1:numel(resultEnv)
            allKapp = [allKapp; resultEnv(i).KappHist(:)]; %#ok<AGROW>
            allKtip = [allKtip; resultEnv(i).KtipHist(:)]; %#ok<AGROW>
        end

        validAll = isfinite(allKapp) & isfinite(allKtip);
        allKapp = allKapp(validAll);
        allKtip = allKtip(validAll);

        xMin = min(allKapp);
        xMax = max(allKapp);
        yMin = min(allKtip);
        yMax = max(allKtip);

        dx = xMax - xMin;
        dy = yMax - yMin;
        if dx <= 0, dx = 1; end
        if dy <= 0, dy = 1; end
        xLimEnv = [xMin - 0.05*dx, xMax + 0.05*dx];
        yLimEnv = [yMin - 0.05*dy, yMax + 0.05*dy];

        nCasesEnv = numel(resultEnv);
        nRowsEnv = ceil(nCasesEnv / nColsEnvelope);

        fig3 = figure('Color','w','Position',[140 70 1200 450*nRowsEnv]);
        tl = tiledlayout(nRowsEnv, nColsEnvelope, 'TileSpacing','compact', 'Padding','compact');

        for i = 1:nCasesEnv
            ax = nexttile(tl); hold(ax,'on'); box(ax,'on');

            if showRawScatterInEnvelope
                plot(ax, resultEnv(i).KappHist, resultEnv(i).KtipHist, '.', ...
                    'Color', [0.72 0.72 0.72], ...
                    'MarkerSize', 6, ...
                    'HandleVisibility','off');
            end

            plot(ax, resultEnv(i).envX, resultEnv(i).envY, '-', ...
                'Color', [0 0.4470 0.7410], ...
                'LineWidth', 2.0, ...
                'DisplayName', 'Upper envelope');

            xlabel(ax, 'K_{app} (MPa\cdotm^{0.5})');
            ylabel(ax, 'K_{tip} (MPa\cdotm^{0.5})');
            title(ax, sprintf('%s  (LL = %.4g)', resultEnv(i).caseName, resultEnv(i).LL), ...
                'Interpreter', 'none');

            xlim(ax, xLimEnv);
            ylim(ax, yLimEnv);

            set(ax, 'FontSize', 12, 'LineWidth', 1.1);
        end

        title(tl, 'K_{app}-K_{tip} upper envelopes for all cases', ...
            'FontSize', 15, 'FontWeight', 'bold');

        if saveEnvelopeData
            envDir = fullfile(figDir, 'envelope_data');
            if ~exist(envDir, 'dir')
                mkdir(envDir);
            end

            for i = 1:nCasesEnv
                envTable = table(resultEnv(i).envX(:), resultEnv(i).envY(:), ...
                    'VariableNames', {'Kapp_MPam05','Ktip_upperEnvelope_MPam05'});
                outEnvCsv = fullfile(envDir, sprintf('%s_Kapp_Ktip_envelope.csv', resultEnv(i).caseName));
                writetable(envTable, outEnvCsv);
            end
            fprintf('包络线数据已保存到: %s\n', envDir);
        end
    end
end

%% ================= 保存图片 =================
if saveFigures
    exportgraphics(fig1, fullfile(figDir, 'fracture_Kapp_N_vs_LL_linear.png'), 'Resolution', 300);
    exportgraphics(fig2, fullfile(figDir, 'fracture_Kapp_N_vs_LL_loglog.png'), 'Resolution', 300);

    if makeEnvelopeFigure && exist('fig3','var') && isgraphics(fig3)
        exportgraphics(fig3, fullfile(figDir, 'Kapp_Ktip_upper_envelopes_subplots.png'), 'Resolution', 300);
    end

    fprintf('图片已保存到: %s\n', figDir);
end

%% ================= 局部函数 =================
function LL = parseLLfromFolderName(folderName)
% 支持:
%   sizeLL2e3  -> LL = 2e3
%   sizeLL8e3  -> LL = 8e3
%   size32e3   -> LL = 32e3

    LL = NaN;

    tok = regexp(folderName, '^sizeLL(\d+)e3$', 'tokens', 'once');
    if ~isempty(tok)
        LL = str2double(tok{1}) * 1e3;
        return;
    end

    tok = regexp(folderName, '^size(\d+)e3$', 'tokens', 'once');
    if ~isempty(tok)
        LL = str2double(tok{1}) * 1e3;
        return;
    end
end

function ktipColName = findKtipColumn(varNames)
% 兼容多种可能的 Ktip 列名
    candidates = { ...
        'Ktip_MPam05', ...
        'Ktip', ...
        'K_tip_MPam05', ...
        'K_tip', ...
        'Ktip_MPam0_5', ...
        'Ktip_MPa_m05'};

    ktipColName = '';
    for i = 1:numel(candidates)
        if ismember(candidates{i}, varNames)
            ktipColName = candidates{i};
            return;
        end
    end
end

function [envX, envY] = computeUpperEnvelope(x, y, method, nBins, minPts, doSmooth, smoothWindow)
% 计算 Kapp-Ktip 上包络线
% 当前实现：
%   method = 'binmax' : 按 x 分箱，每箱取最大 y

    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    envX = [];
    envY = [];

    if numel(x) < 2
        return;
    end

    switch lower(method)
        case 'binmax'
            xMin = min(x);
            xMax = max(x);
            if xMax <= xMin
                envX = x;
                envY = y;
                return;
            end

            edges = linspace(xMin, xMax, nBins+1);

            for i = 1:nBins
                if i < nBins
                    idx = x >= edges(i) & x < edges(i+1);
                else
                    idx = x >= edges(i) & x <= edges(i+1);
                end

                if nnz(idx) < minPts
                    continue;
                end

                xi = x(idx);
                yi = y(idx);

                [ymax, imax] = max(yi);
                xmax = xi(imax);

                envX(end+1,1) = xmax; %#ok<AGROW>
                envY(end+1,1) = ymax; %#ok<AGROW>
            end

            if isempty(envX)
                return;
            end

            [envX, ord] = sort(envX);
            envY = envY(ord);

            % 去除重复 x
            [envX, ia] = unique(envX, 'stable');
            envY = envY(ia);

            if doSmooth && numel(envY) >= smoothWindow
                envY = movmean(envY, smoothWindow);
            end

        otherwise
            error('未知包络线方法: %s', method);
    end
end