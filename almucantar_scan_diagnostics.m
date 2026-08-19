%% almucantar_scan_diagnostics.m
% RAW AERONET almucantar branch-symmetry diagnostics.
%
% IMPORTANT - SDAT CHECK:
% This script is NOT intended to reproduce the photometer subset written
% to the GRASP SDAT file. It analyses the full raw AERONET almucantar scan,
% including scattering-angle information and both branches.
%
% For a direct check against Mariana's SDAT construction, use:
%
%     almucantar_sdat_crosscheck.m
%
% That companion script reproduces the 30 positive/outward azimuths,
% applies the common four-wavelength NaN/missing-value mask, and reports
% N_SDAT_common, retained azimuths, RAA, SZA and VZA. It deliberately does
% not use scattering angle because scattering angle is not read when the
% SDAT is constructed.
%
% This file is retained for the separate scientific diagnostic of whether
% the raw almucantar scan shows branch asymmetry / possible horizontal
% inhomogeneity. The 20% line is a reference diagnostic only, not a GRASP
% QC threshold.
%
% Branch relative difference:
%
%     100 * abs(I_plus - I_minus) / ((I_plus + I_minus)/2)
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

% Reference line/flag only; not a proposed selection threshold.
branchDifferenceThresholdPct = 20;

% Default target cases for the July/August diagnostic file.
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

    complete = all(ismember(graspWavelengths, wls));

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
                     targetToleranceMinutes, string(targetTimes(it)), dmin);
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

        valid = isfinite(rad) & isfinite(sca) & rad > -900 & sca > -900;
        absNom = abs(nominalAngles);

        Ncrit_2to6   = sum(valid & absNom >= 2  & absNom <= 6);
        Ncrit_6to30  = sum(valid & absNom >= 6  & absNom <= 30);
        Ncrit_30to80 = sum(valid & absNom >= 30 & absNom <= 80);
        Ncrit_gt80   = sum(valid & absNom > 80);

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

        [pairNom, pairScat, Iplus, Iminus, branchDiff] = ...
            branchSymmetry(nominalAngles, rad, sca, branchMinScatteringAngle);

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
        waveRows(iw).Ncrit_2to6 = Ncrit_2to6;
        waveRows(iw).Ncrit_6to30 = Ncrit_6to30;
        waveRows(iw).Ncrit_30to80 = Ncrit_30to80;
        waveRows(iw).Ncrit_gt80 = Ncrit_gt80;
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
    inputFile = char(inputFile);
    cleanupObj = [];

    if endsWith(lower(inputFile), '.zip')
        tmpDir = tempname;
        mkdir(tmpDir);
        unzip(inputFile, tmpDir);
        d = dir(fullfile(tmpDir, '**', '*.alm'));
        if isempty(d)
            error('ZIP file does not contain an .alm file.');
        elseif numel(d) > 1
            error('ZIP contains more than one .alm file; extract and specify one file.');
        end
        almFile = fullfile(d(1).folder, d(1).name);
        cleanupObj = onCleanup(@() rmdir(tmpDir,'s'));
    elseif endsWith(lower(inputFile), '.alm')
        almFile = inputFile;
    else
        error('inputFile must be an .alm file or .zip containing one .alm file.');
    end
end

function [pairNom, pairScat, Iplus, Iminus, branchDiff] = ...
    branchSymmetry(nominalAngles, rad, sca, minScat)

    tol = 1e-8;
    positiveAngles = unique(nominalAngles(nominalAngles > 0));

    pairNom = [];
    pairScat = [];
    Iplus = [];
    Iminus = [];
    branchDiff = [];

    for ia = 1:numel(positiveAngles)
        a = positiveAngles(ia);
        ip = find(abs(nominalAngles-a) < tol);
        im = find(abs(nominalAngles+a) < tol);

        if isempty(ip) || isempty(im)
            continue
        end

        % If an angle occurs more than once in the raw scan, pair the first
        % valid occurrence on each branch for this diagnostic.
        ip = ip(isfinite(rad(ip)) & isfinite(sca(ip)) & rad(ip) > -900 & sca(ip) > -900);
        im = im(isfinite(rad(im)) & isfinite(sca(im)) & rad(im) > -900 & sca(im) > -900);
        if isempty(ip) || isempty(im)
            continue
        end
        ip = ip(1);
        im = im(1);

        scatPair = mean(abs([sca(ip),sca(im)]));
        if scatPair <= minScat
            continue
        end

        p = rad(ip);
        m = rad(im);
        denom = (p+m)/2;
        if ~isfinite(denom) || denom <= 0
            continue
        end

        pairNom(end+1,1) = a; %#ok<AGROW>
        pairScat(end+1,1) = scatPair; %#ok<AGROW>
        Iplus(end+1,1) = p; %#ok<AGROW>
        Iminus(end+1,1) = m; %#ok<AGROW>
        branchDiff(end+1,1) = 100*abs(p-m)/denom; %#ok<AGROW>
    end
end
