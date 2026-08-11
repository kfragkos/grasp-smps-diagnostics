%% almucantar_scan_diagnostics.m
% Diagnostics for AERONET Version 3 RAW ALMUCANTAR (.alm) files.
%
% The script reproduces the diagnostics discussed for the GRASP-SMPS
% comparison:
%   1) number of valid raw angular radiance measurements at 440, 675,
%      870 and 1020 nm;
%   2) counts in the manuscript's nominal azimuth-angle intervals
%      [2,6], [6,30], [30,80], >80 deg (these overlap at 6 and 30 deg,
%      so they are NOT additive), plus a non-overlapping partition for
%      bookkeeping;
%   3) solar zenith angle (SZA);
%   4) maximum ACTUAL scattering angle represented in the raw scan;
%   5) symmetry between the two almucantar branches for corresponding
%      +/- nominal angles, restricted to actual scattering angle > 6 deg;
%   6) mean, median and maximum branch difference and the number of pairs;
%   7) scan-level summaries across the four GRASP wavelengths.
%
% IMPORTANT:
% This analyses the RAW .alm measurements. "NValid" therefore means
% valid points present in the raw AERONET scan. If additional points were
% removed before writing the GRASP SDAT file, compare against the SDAT or
% the exact preprocessing code before calling them the points "supplied
% to GRASP".
%
% Branch relative difference used here:
%
%     100 * abs(I_plus - I_minus) / ((I_plus + I_minus)/2)
%
% This is the same definition used to reproduce the July/August numbers
% discussed in the manuscript analysis.
%
% No special MATLAB toolboxes are required.

clear; clc;

%% ------------------------ USER SETTINGS -------------------------------

% Can be either the .zip downloaded from AERONET or the extracted .alm file.
inputFile = '20240710_20240815_Magurele_Inoe.zip';

% GRASP wavelengths used in this study.
graspWavelengths = [440 675 870 1020];

% Rows from one almucantar scan are separated by ~1-2 min, whereas
% consecutive scans are usually separated by much longer. A gap >5 min
% starts a new scan.
scanGapMinutes = 5;

% Maximum allowed difference between a requested target time and the
% nearest complete four-wavelength scan.
targetToleranceMinutes = 8;

% Branch-symmetry diagnostics are calculated only for corresponding branch
% pairs with actual scattering angle > 6 deg.
branchMinScatteringAngle = 6;

% Useful reference line/flag. Change if you want to test another value.
branchDifferenceThresholdPct = 20;

% Target cases to inspect.
% Times are approximate scan-centre times; the code finds the nearest
% complete 440/675/870/1020-nm scan.
targetTimes = [
    datetime(2024,7,13,14,35,0,'TimeZone','UTC')
    datetime(2024,7,13,14,54,0,'TimeZone','UTC')
    datetime(2024,8,12,14,11,0,'TimeZone','UTC')
    datetime(2024,8,12,14,46,0,'TimeZone','UTC')
    ];

targetLabels = [
    "13Jul_good"
    "13Jul_poor"
    "12Aug_poor"
    "12Aug_good"
    ];

% Set to [] to analyse ALL complete scans instead of selected target times:
% targetTimes = [];
% targetLabels = strings(0,1);

makePlots = true;

% Output folder. Empty -> folder containing input file.
outDir = fileparts(inputFile);
if isempty(outDir)
    outDir = pwd;
end

%% -------------------------- READ FILE --------------------------------

[almFile, cleanupObj] = getAlmFile(inputFile); %#ok<NASGU>

lines = readlines(almFile);

headerIdx = find(startsWith(lines, "AERONET_Site,Date("), 1, 'first');
if isempty(headerIdx)
    error('Could not find the AERONET data header in %s', almFile);
end

headerFields = split(lines(headerIdx), ',');
if strlength(headerFields(end)) == 0
    headerFields(end) = [];
end

nMeta = 11;
nRemaining = numel(headerFields) - nMeta;

if mod(nRemaining, 2) ~= 0
    error(['Unexpected .alm format: after the first %d metadata columns, ' ...
           'the remaining columns cannot be divided equally into radiance ' ...
           'and scattering-angle fields.'], nMeta);
end

nAngles = nRemaining / 2;

nominalAngles = str2double(headerFields(nMeta+1:nMeta+nAngles));
nominalAngles = nominalAngles(:);

% Verify that the second half of the numeric header repeats the nominal
% angle sequence for the scattering-angle fields.
nominalAngles2 = str2double(headerFields(nMeta+nAngles+1:nMeta+2*nAngles));
nominalAngles2 = nominalAngles2(:);

if any(abs(nominalAngles - nominalAngles2) > 1e-10)
    warning('The radiance/scattering-angle header sequences are not identical.');
end

% Parse only wavelengths used by GRASP.
rec = struct('Time',{},'Wavelength_nm',{},'SZA_deg',{}, ...
             'Radiance',{},'ScatteringAngle_deg',{});

k = 0;
for iline = headerIdx+1:numel(lines)

    thisLine = strtrim(lines(iline));
    if strlength(thisLine) == 0
        continue
    end

    f = split(thisLine, ',');

    % Remove possible final blank field caused by a trailing comma.
    if strlength(f(end)) == 0
        f(end) = [];
    end

    if numel(f) < nMeta + 2*nAngles
        continue
    end

    wl = str2double(f(5));
    if ~ismember(round(wl), graspWavelengths)
        continue
    end

    try
        t = datetime(f(2) + " " + f(3), ...
            'InputFormat','dd:MM:yyyy HH:mm:ss', ...
            'TimeZone','UTC');
    catch
        continue
    end

    sza = str2double(f(7));

    rad = str2double(f(nMeta+1:nMeta+nAngles));
    sca = str2double(f(nMeta+nAngles+1:nMeta+2*nAngles));

    k = k + 1;
    rec(k).Time = t;
    rec(k).Wavelength_nm = round(wl);
    rec(k).SZA_deg = sza;
    rec(k).Radiance = rad(:);
    rec(k).ScatteringAngle_deg = sca(:);
end

if isempty(rec)
    error('No 440/675/870/1020-nm almucantar rows were found.');
end

% Sort rows chronologically.
[~,ord] = sort([rec.Time]);
rec = rec(ord);

%% -------------------------- BUILD SCANS -------------------------------

% A new scan begins when the time gap between successive selected
% wavelength rows exceeds scanGapMinutes.
nRec = numel(rec);
scanID = ones(nRec,1);

for i = 2:nRec
    dtmin = minutes(rec(i).Time - rec(i-1).Time);
    if dtmin > scanGapMinutes
        scanID(i) = scanID(i-1) + 1;
    else
        scanID(i) = scanID(i-1);
    end
end

allScanIDs = unique(scanID);

scan = struct('ID',{},'CenterTime',{},'RecordIndices',{},'Complete',{});
ks = 0;

for is = 1:numel(allScanIDs)
    idx = find(scanID == allScanIDs(is));
    wls = [rec(idx).Wavelength_nm];

    % Require at least one row for all four GRASP wavelengths.
    complete = all(ismember(graspWavelengths, wls));

    % Use only the GRASP wavelengths for scan centre.
    idx4 = idx(ismember(wls, graspWavelengths));
    tt = [rec(idx4).Time];
    tsec = posixtime(tt);
    centerTime = datetime(mean(tsec), 'ConvertFrom','posixtime', ...
                          'TimeZone','UTC');

    ks = ks + 1;
    scan(ks).ID = allScanIDs(is);
    scan(ks).CenterTime = centerTime;
    scan(ks).RecordIndices = idx;
    scan(ks).Complete = complete;
end

completeIdx = find([scan.Complete]);

if isempty(completeIdx)
    error('No complete 440/675/870/1020-nm scans were found.');
end

%% ------------------------- SELECT SCANS -------------------------------

selectedScanIdx = [];
selectedTargetIndex = [];
selectedTargetTime = datetime.empty(0,1);
selectedLabel = strings(0,1);

if isempty(targetTimes)
    selectedScanIdx = completeIdx(:);
    selectedTargetIndex = (1:numel(selectedScanIdx))';
    selectedTargetTime = reshape([scan(selectedScanIdx).CenterTime],[],1);
    selectedLabel = "scan_" + string(selectedTargetIndex);
else
    if numel(targetLabels) ~= numel(targetTimes)
        error('targetLabels and targetTimes must have the same number of entries.');
    end

    centers = reshape([scan(completeIdx).CenterTime],[],1);

    for it = 1:numel(targetTimes)
        deltaMin = abs(minutes(centers - targetTimes(it)));
        [dmin,j] = min(deltaMin);

        if dmin > targetToleranceMinutes
            warning(['No complete scan was found within %.1f min of %s. ' ...
                     'Nearest scan is %.2f min away.'], ...
                     targetToleranceMinutes, ...
                     string(targetTimes(it)), dmin);
            continue
        end

        selectedScanIdx(end+1,1) = completeIdx(j); %#ok<SAGROW>
        selectedTargetIndex(end+1,1) = it; %#ok<SAGROW>
        selectedTargetTime(end+1,1) = targetTimes(it); %#ok<SAGROW>
        selectedLabel(end+1,1) = targetLabels(it); %#ok<SAGROW>
    end
end

if isempty(selectedScanIdx)
    error('No requested scans could be matched.');
end

%% ----------------------- CALCULATE DIAGNOSTICS ------------------------

waveRows = struct([]);
pairRows = struct([]);

iw = 0;
ipair = 0;

for isel = 1:numel(selectedScanIdx)

    s = scan(selectedScanIdx(isel));
    idx = s.RecordIndices;

    % Keep exactly one record at each requested wavelength.
    % If duplicates ever occur, take the one closest to scan centre.
    thisWlMeans = nan(numel(graspWavelengths),1);
    thisWlMax = nan(numel(graspWavelengths),1);

    for jw = 1:numel(graspWavelengths)

        wlWanted = graspWavelengths(jw);
        candidates = idx([rec(idx).Wavelength_nm] == wlWanted);

        if isempty(candidates)
            continue
        end

        if numel(candidates) > 1
            d = abs(minutes([rec(candidates).Time] - s.CenterTime));
            [~,jj] = min(d);
            ir = candidates(jj);
        else
            ir = candidates(1);
        end

        rad = rec(ir).Radiance;
        sca = rec(ir).ScatteringAngle_deg;

        % AERONET missing values are typically -999 or -999.000000.
        valid = isfinite(rad) & isfinite(sca) & rad > -900 & sca > -900;

        absNom = abs(nominalAngles);

        % Counts in the exact intervals stated in the manuscript.
        % NOTE: [2,6], [6,30] and [30,80] overlap at the boundaries,
        % therefore these four counts must NOT be summed to recover NValid.
        Ncrit_2to6   = sum(valid & absNom >= 2  & absNom <= 6);
        Ncrit_6to30  = sum(valid & absNom >= 6  & absNom <= 30);
        Ncrit_30to80 = sum(valid & absNom >= 30 & absNom <= 80);
        Ncrit_gt80   = sum(valid & absNom > 80);

        % Non-overlapping partition of the same valid measurements.
        % These four DO sum to NValid.
        Npart_2to6   = sum(valid & absNom >= 2  & absNom <= 6);
        Npart_6to30  = sum(valid & absNom >  6  & absNom <= 30);
        Npart_30to80 = sum(valid & absNom > 30  & absNom <= 80);
        Npart_gt80   = sum(valid & absNom > 80);

        NValid = sum(valid);

        if any(valid)
            maxActualScat = max(abs(sca(valid)));
        else
            maxActualScat = NaN;
        end

        % Pair corresponding +angle and -angle measurements and calculate
        % branch symmetry for actual scattering angle > 6 deg.
        [pairNom, pairScat, Iplus, Iminus, branchDiff] = ...
            branchSymmetry(nominalAngles, rad, sca, ...
                           branchMinScatteringAngle);

        if isempty(branchDiff)
            nPairs = 0;
            meanBranch = NaN;
            medianBranch = NaN;
            maxBranch = NaN;
            nAboveThreshold = 0;
        else
            nPairs = numel(branchDiff);
            meanBranch = mean(branchDiff);
            medianBranch = median(branchDiff);
            maxBranch = max(branchDiff);
            nAboveThreshold = sum(branchDiff > branchDifferenceThresholdPct);
        end

        iw = iw + 1;
        waveRows(iw).TargetIndex = selectedTargetIndex(isel);
        waveRows(iw).Label = selectedLabel(isel);
        waveRows(iw).RequestedTime_UTC = selectedTargetTime(isel);
        waveRows(iw).ScanCenter_UTC = s.CenterTime;
        waveRows(iw).RowTime_UTC = rec(ir).Time;
        waveRows(iw).Wavelength_nm = wlWanted;
        waveRows(iw).SZA_deg = rec(ir).SZA_deg;
        waveRows(iw).NValidRaw = NValid;

        % Exact manuscript intervals (overlapping boundaries; not additive).
        waveRows(iw).Ncrit_2to6 = Ncrit_2to6;
        waveRows(iw).Ncrit_6to30 = Ncrit_6to30;
        waveRows(iw).Ncrit_30to80 = Ncrit_30to80;
        waveRows(iw).Ncrit_gt80 = Ncrit_gt80;

        % Non-overlapping partition (additive to NValidRaw).
        waveRows(iw).Npart_2to6 = Npart_2to6;
        waveRows(iw).Npart_6to30 = Npart_6to30;
        waveRows(iw).Npart_30to80 = Npart_30to80;
        waveRows(iw).Npart_gt80 = Npart_gt80;
        waveRows(iw).MaxActualScatteringAngle_deg = maxActualScat;
        waveRows(iw).NBranchPairs_gt6deg = nPairs;
        waveRows(iw).MeanBranchDiff_pct = meanBranch;
        waveRows(iw).MedianBranchDiff_pct = medianBranch;
        waveRows(iw).MaxBranchDiff_pct = maxBranch;
        waveRows(iw).NPairsAbove20pct = nAboveThreshold;
        waveRows(iw).BranchMaxBelow20pct = ...
            (~isnan(maxBranch) && maxBranch <= branchDifferenceThresholdPct);

        thisWlMeans(jw) = meanBranch;
        thisWlMax(jw) = maxBranch;

        % Save each individual branch pair.
        for jp = 1:numel(branchDiff)
            ipair = ipair + 1;
            pairRows(ipair).TargetIndex = selectedTargetIndex(isel);
            pairRows(ipair).Label = selectedLabel(isel);
            pairRows(ipair).ScanCenter_UTC = s.CenterTime;
            pairRows(ipair).RowTime_UTC = rec(ir).Time;
            pairRows(ipair).Wavelength_nm = wlWanted;
            pairRows(ipair).NominalAbsAngle_deg = pairNom(jp);
            pairRows(ipair).ActualScatteringAngle_deg = pairScat(jp);
            pairRows(ipair).RadiancePlus = Iplus(jp);
            pairRows(ipair).RadianceMinus = Iminus(jp);
            pairRows(ipair).BranchDiff_pct = branchDiff(jp);
        end
    end
end

waveTable = struct2table(waveRows);
pairTable = struct2table(pairRows);

% Scan-level summary generated from wavelength table.
scanRows = struct([]);

for isel = 1:numel(selectedScanIdx)

    mask = waveTable.TargetIndex == selectedTargetIndex(isel);

    scanRows(isel).TargetIndex = selectedTargetIndex(isel);
    scanRows(isel).Label = selectedLabel(isel);
    scanRows(isel).RequestedTime_UTC = selectedTargetTime(isel);
    scanRows(isel).ScanCenter_UTC = scan(selectedScanIdx(isel)).CenterTime;
    scanRows(isel).MeanSZA_deg = mean(waveTable.SZA_deg(mask),'omitnan');
    scanRows(isel).MinNValidRawAcross4WL = min(waveTable.NValidRaw(mask));
    scanRows(isel).MaxActualScatteringAngle_deg = ...
        max(waveTable.MaxActualScatteringAngle_deg(mask),[],'omitnan');

    % This definition reproduces the scan-average branch statistics used
    % previously: mean of the four wavelength-specific means.
    scanRows(isel).MeanOfWavelengthMeanBranchDiff_pct = ...
        mean(waveTable.MeanBranchDiff_pct(mask),'omitnan');

    scanRows(isel).MaximumBranchDiffAcross4WL_pct = ...
        max(waveTable.MaxBranchDiff_pct(mask),[],'omitnan');

    scanRows(isel).All4WavelengthsMaxBelow20pct = ...
        all(waveTable.MaxBranchDiff_pct(mask) <= branchDifferenceThresholdPct);
end

scanTable = struct2table(scanRows);

%% --------------------------- DISPLAY ---------------------------------

fprintf('\n================ SCAN-LEVEL SUMMARY ================\n');
disp(scanTable);

fprintf('\n============= WAVELENGTH-LEVEL DIAGNOSTICS =============\n');
disp(waveTable(:, { ...
    'Label','ScanCenter_UTC','RowTime_UTC','Wavelength_nm','SZA_deg', ...
    'NValidRaw','Ncrit_2to6','Ncrit_6to30','Ncrit_30to80','Ncrit_gt80', ...
    'MaxActualScatteringAngle_deg','NBranchPairs_gt6deg', ...
    'MeanBranchDiff_pct','MedianBranchDiff_pct','MaxBranchDiff_pct', ...
    'NPairsAbove20pct'}));

%% -------------------------- SAVE TABLES -------------------------------

waveCsv = fullfile(outDir, 'almucantar_wavelength_diagnostics.csv');
scanCsv = fullfile(outDir, 'almucantar_scan_summary.csv');
pairCsv = fullfile(outDir, 'almucantar_branch_pairs.csv');

writetable(waveTable, waveCsv);
writetable(scanTable, scanCsv);
writetable(pairTable, pairCsv);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n', waveCsv, scanCsv, pairCsv);

%% -------------------------- OPTIONAL PLOTS ----------------------------

if makePlots

    plotWavelengths = [440 675 870 1020];

    figure('Color','w','Name','Almucantar branch-symmetry diagnostics');

    for jw = 1:numel(plotWavelengths)
        subplot(2,2,jw);
        hold on; box on;

        wl = plotWavelengths(jw);

        for isel = 1:numel(selectedScanIdx)
            m = pairTable.TargetIndex == selectedTargetIndex(isel) & ...
                pairTable.Wavelength_nm == wl;

            if any(m)
                [x,oo] = sort(pairTable.ActualScatteringAngle_deg(m));
                yy = pairTable.BranchDiff_pct(m);
                yy = yy(oo);

                plot(x, yy, '-o', 'LineWidth', 1.0, 'MarkerSize', 4, ...
                     'DisplayName', char(selectedLabel(isel)));
            end
        end

        yline(branchDifferenceThresholdPct,'k--','20%');
        xlabel('Actual scattering angle [deg]');
        ylabel('Branch difference [%]');
        title(sprintf('%d nm',wl));
        legend('Location','best');
        grid on;
    end

    figFile = fullfile(outDir,'almucantar_branch_symmetry.png');
    exportgraphics(gcf, figFile, 'Resolution', 180);
    fprintf('  %s\n', figFile);
end

%% ========================== LOCAL FUNCTIONS ===========================

function [almFile, cleanupObj] = getAlmFile(inputFile)
% Accept either a ZIP containing one .alm file or an extracted .alm file.

    inputFile = char(inputFile);
    cleanupObj = [];

    if endsWith(lower(inputFile), '.zip')
        tmpDir = tempname;
        mkdir(tmpDir);

        files = unzip(inputFile, tmpDir);
        isAlm = endsWith(lower(string(files)), '.alm');

        if ~any(isAlm)
            rmdir(tmpDir,'s');
            error('ZIP file contains no .alm file.');
        end

        almFiles = string(files(isAlm));
        if numel(almFiles) > 1
            warning('More than one .alm file found; using the first.');
        end

        almFile = char(almFiles(1));
        cleanupObj = onCleanup(@() safeRemoveDir(tmpDir));

    elseif endsWith(lower(inputFile), '.alm')
        almFile = inputFile;

        if ~isfile(almFile)
            error('Input .alm file does not exist: %s',almFile);
        end
    else
        error('inputFile must be either an AERONET .alm file or a .zip containing one.');
    end
end

function safeRemoveDir(d)
    if isfolder(d)
        rmdir(d,'s');
    end
end

function [nomAbs, scatPair, Iplus, Iminus, diffPct] = ...
    branchSymmetry(nominalAngles, rad, scat, minScat)
% Pair the two outer almucantar branches at equal absolute nominal angles.
%
% Only |nominal angle| > 6 deg is considered. A pair is retained if the
% average absolute ACTUAL scattering angle of the +/- measurements exceeds
% minScat.
%
% The returned branch difference is:
% 100*abs(Iplus-Iminus)/mean([Iplus,Iminus]).

    nominalAngles = nominalAngles(:);
    rad = rad(:);
    scat = scat(:);

    valid = isfinite(rad) & isfinite(scat) & rad > -900 & scat > -900;

    outerAngles = unique(abs(nominalAngles(abs(nominalAngles) > 6)));
    outerAngles = sort(outerAngles);

    nomAbs = [];
    scatPair = [];
    Iplus = [];
    Iminus = [];
    diffPct = [];

    for j = 1:numel(outerAngles)

        a = outerAngles(j);

        ip = find(nominalAngles ==  a, 1, 'first');
        im = find(nominalAngles == -a, 1, 'first');

        if isempty(ip) || isempty(im)
            continue
        end

        if ~(valid(ip) && valid(im))
            continue
        end

        thisScat = mean(abs([scat(ip), scat(im)]));

        if thisScat <= minScat
            continue
        end

        denom = (rad(ip) + rad(im)) / 2;
        if denom <= 0
            continue
        end

        thisDiff = 100 * abs(rad(ip) - rad(im)) / denom;

        nomAbs(end+1,1) = a; %#ok<AGROW>
        scatPair(end+1,1) = thisScat; %#ok<AGROW>
        Iplus(end+1,1) = rad(ip); %#ok<AGROW>
        Iminus(end+1,1) = rad(im); %#ok<AGROW>
        diffPct(end+1,1) = thisDiff; %#ok<AGROW>
    end
end
