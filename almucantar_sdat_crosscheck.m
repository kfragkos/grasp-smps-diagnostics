%% almucantar_sdat_crosscheck.m
% Cross-check the AERONET almucantar radiances retained when constructing
% a GRASP SDAT file.
%
% The script reproduces the photometer-side angular selection:
%   1) read an annual AERONET raw almucantar .alm file;
%   2) select scans containing all four GRASP wavelengths
%      (440, 675, 870, 1020 nm);
%   3) use the 30 positive/outward almucantar azimuth positions from
%      2 deg to 180 deg (the repeated 6 deg entry is retained once);
%   4) report the number of valid radiances at each wavelength and the
%      angular-range coverage;
%   5) apply a common four-wavelength mask: if a radiance is missing/NaN
%      at any wavelength, that azimuth is removed from all four;
%   6) report retained azimuths, RAA = azimuth + 180, SZA and
%      VZA = 180 - SZA for comparison with an independently generated SDAT.
%
% IMPORTANT:
% - This cross-check deliberately does NOT use the scattering-angle fields
%   from the .alm file because they are not required for this SDAT step.
% - Raw radiances are included only to trace the common mask. They are not
%   normalized here by E0/Earth-Sun distance.
% - Branch-symmetry/scattering-angle diagnostics are handled separately in
%   almucantar_scan_diagnostics.m.
%
% No special MATLAB toolboxes are required.

clear; clc;

%% ------------------------ USER SETTINGS -------------------------------

% Annual AERONET raw almucantar file. A ZIP containing one .alm file is
% also accepted.
inputFile = '20240101_20241231_Magurele_Inoe.alm';

% GRASP photometer wavelengths used in the SDAT.
graspWavelengths = [440 675 870 1020];

% Consecutive wavelength rows belonging to one almucantar scan are normally
% separated by ~1-2 min. A gap >5 min starts a new scan.
scanGapMinutes = 5;

% Maximum difference between a requested retrieval time and the nearest
% complete four-wavelength almucantar scan.
targetToleranceMinutes = 8;

% Example retrieval times from the study. Edit as needed.
targetTimes = [
    datetime(2024,5,26,14,26,03,'TimeZone','UTC')
    datetime(2024,5,26,16,00,28,'TimeZone','UTC')
    datetime(2024,5,26,16,19,27,'TimeZone','UTC')
    datetime(2024,7,13,14,35,09,'TimeZone','UTC')
    datetime(2024,7,13,14,53,53,'TimeZone','UTC')
    datetime(2024,7,14,14,35,15,'TimeZone','UTC')
    datetime(2024,7,15,14,54,02,'TimeZone','UTC')
    datetime(2024,8,08,16,08,24,'TimeZone','UTC')
    datetime(2024,8,11,14,12,35,'TimeZone','UTC')
    datetime(2024,8,11,14,46,50,'TimeZone','UTC')
    datetime(2024,9,04,13,47,49,'TimeZone','UTC')
    datetime(2024,9,23,07,19,17,'TimeZone','UTC')
    ];

targetLabels = string(datestr(targetTimes,'yyyymmdd_HHMMSS'));

% To analyse every complete scan instead of selected retrieval times:
% targetTimes = [];
% targetLabels = strings(0,1);

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
if mod(nRemaining,2) ~= 0
    error('Unexpected .alm structure: radiance/scattering fields cannot be split equally.');
end

% The raw file stores one block of radiances and one block of scattering
% angles. Only the radiance block is used here.
nAnglesRaw = nRemaining/2;
nominalAnglesRaw = str2double(headerFields(nMeta+1:nMeta+nAnglesRaw));
nominalAnglesRaw = nominalAnglesRaw(:);

%% Identify the 30 positive/outward azimuth positions used for SDAT

iFirst2 = find(abs(nominalAnglesRaw-2) < 1e-10, 1, 'first');
iFirst180 = find(abs(nominalAnglesRaw-180) < 1e-10, 1, 'first');

if isempty(iFirst2) || isempty(iFirst180) || iFirst180 <= iFirst2
    error('Could not identify the positive outward +2 to +180 deg branch in the .alm header.');
end

candidateIdx = (iFirst2:iFirst180)';
candidateAngles = nominalAnglesRaw(candidateIdx);

% Stable de-duplication of the positive azimuth sequence.
keepUnique = false(size(candidateAngles));
seenAngles = [];
for i = 1:numel(candidateAngles)
    if candidateAngles(i) > 0 && ~any(abs(seenAngles-candidateAngles(i)) < 1e-10)
        keepUnique(i) = true;
        seenAngles(end+1,1) = candidateAngles(i); %#ok<SAGROW>
    end
end

sdatRawIdx = candidateIdx(keepUnique);
sdatAzimuth_deg = nominalAnglesRaw(sdatRawIdx);  % column vector
sdatAzRow = sdatAzimuth_deg';                    % row vector

if numel(sdatRawIdx) ~= 30
    warning('Expected 30 initial SDAT azimuths, but identified %d.', numel(sdatRawIdx));
end

fprintf('Initial SDAT azimuth positions identified: %d\n', numel(sdatRawIdx));
fprintf('%s\n', join(string(sdatAzRow), ', '));

%% Parse the four photometer wavelengths; do NOT read scattering angle

rec = struct('Time',{},'Wavelength_nm',{},'SZA_deg',{},'Radiance',{});
k = 0;

for iline = headerIdx+1:numel(lines)
    thisLine = strtrim(lines(iline));
    if strlength(thisLine) == 0
        continue
    end

    f = split(thisLine, ',');
    if strlength(f(end)) == 0
        f(end) = [];
    end
    if numel(f) < nMeta + nAnglesRaw
        continue
    end

    wl = str2double(f(5));
    if ~ismember(round(wl), graspWavelengths)
        continue
    end

    try
        t = datetime(f(2) + " " + f(3), ...
            'InputFormat','dd:MM:yyyy HH:mm:ss', 'TimeZone','UTC');
    catch
        continue
    end

    k = k + 1;
    rec(k).Time = t;
    rec(k).Wavelength_nm = round(wl);
    rec(k).SZA_deg = str2double(f(7));
    rec(k).Radiance = str2double(f(nMeta+1:nMeta+nAnglesRaw));
    rec(k).Radiance = rec(k).Radiance(:);
end

if isempty(rec)
    error('No 440/675/870/1020-nm almucantar rows found.');
end

[~,ord] = sort([rec.Time]);
rec = rec(ord);

%% -------------------------- BUILD SCANS -------------------------------

nRec = numel(rec);
scanID = ones(nRec,1);
for i = 2:nRec
    if minutes(rec(i).Time-rec(i-1).Time) > scanGapMinutes
        scanID(i) = scanID(i-1)+1;
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
    complete = all(ismember(graspWavelengths,wls));

    tt = [rec(idx).Time];
    centerTime = datetime(mean(posixtime(tt)), ...
        'ConvertFrom','posixtime','TimeZone','UTC');

    ks = ks+1;
    scan(ks).ID = allScanIDs(is);
    scan(ks).CenterTime = centerTime;
    scan(ks).RecordIndices = idx;
    scan(ks).Complete = complete;
end

completeIdx = find([scan.Complete]);
if isempty(completeIdx)
    error('No complete four-wavelength scans found.');
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
        error('targetLabels and targetTimes must contain the same number of entries.');
    end

    centers = reshape([scan(completeIdx).CenterTime],[],1);
    for it = 1:numel(targetTimes)
        [dmin,j] = min(abs(minutes(centers-targetTimes(it))));
        if dmin > targetToleranceMinutes
            warning('No complete scan within %.1f min of %s (nearest %.2f min).', ...
                targetToleranceMinutes,string(targetTimes(it)),dmin);
            continue
        end
        selectedScanIdx(end+1,1) = completeIdx(j); %#ok<SAGROW>
        selectedTargetIndex(end+1,1) = it; %#ok<SAGROW>
        selectedTargetTime(end+1,1) = targetTimes(it); %#ok<SAGROW>
        selectedLabel(end+1,1) = targetLabels(it); %#ok<SAGROW>
    end
end

if isempty(selectedScanIdx)
    error('None of the requested scans could be matched.');
end

%% ----------------------- SDAT CROSS-CHECK -----------------------------

summaryRows = struct([]);
angleRows = struct([]);
ia = 0;

for isel = 1:numel(selectedScanIdx)
    s = scan(selectedScanIdx(isel));
    idxScan = s.RecordIndices;

    rad4 = nan(numel(graspWavelengths),numel(sdatRawIdx));
    sza4 = nan(numel(graspWavelengths),1);
    rowTime4 = NaT(numel(graspWavelengths),1,'TimeZone','UTC');

    % Select exactly one row for each wavelength. If duplicates occur,
    % choose the row nearest the scan centre.
    for jw = 1:numel(graspWavelengths)
        wl = graspWavelengths(jw);
        candidates = idxScan([rec(idxScan).Wavelength_nm] == wl);
        if isempty(candidates)
            continue
        end

        if numel(candidates) > 1
            [~,jj] = min(abs(minutes([rec(candidates).Time]-s.CenterTime)));
            ir = candidates(jj);
        else
            ir = candidates(1);
        end

        rad4(jw,:) = rec(ir).Radiance(sdatRawIdx);
        sza4(jw) = rec(ir).SZA_deg;
        rowTime4(jw) = rec(ir).Time;
    end

    % Validity at each wavelength before the common four-wavelength mask.
    valid4 = isfinite(rad4) & rad4 > -900;
    nValidPerWL = sum(valid4,2);

    % One invalid wavelength removes this azimuth from all four.
    commonKeep = all(valid4,1);  % 1 x N
    retainedAz = sdatAzimuth_deg(commonKeep');
    removedAz = sdatAzimuth_deg(~commonKeep');
    retainedRAA = retainedAz + 180;

    % Non-overlapping angular-range coverage after the common mask.
    n2to6   = sum(commonKeep & sdatAzRow >= 2 & sdatAzRow <= 6);
    n6to30  = sum(commonKeep & sdatAzRow >  6 & sdatAzRow <= 30);
    n30to80 = sum(commonKeep & sdatAzRow > 30 & sdatAzRow <= 80);
    ngt80   = sum(commonKeep & sdatAzRow > 80);

    % Per-wavelength range coverage before the common mask.
    rangeOKperWL = false(4,1);
    for jw = 1:4
        v = valid4(jw,:);
        c1 = sum(v & sdatAzRow >= 2 & sdatAzRow <= 6);
        c2 = sum(v & sdatAzRow >  6 & sdatAzRow <= 30);
        c3 = sum(v & sdatAzRow > 30 & sdatAzRow <= 80);
        c4 = sum(v & sdatAzRow > 80);
        rangeOKperWL(jw) = all([c1 c2 c3 c4] >= 1);
    end

    summaryRows(isel).TargetIndex = selectedTargetIndex(isel);
    summaryRows(isel).Label = selectedLabel(isel);
    summaryRows(isel).RequestedTime_UTC = selectedTargetTime(isel);
    summaryRows(isel).MatchedScanCenter_UTC = s.CenterTime;
    summaryRows(isel).TimeDifference_min = minutes(s.CenterTime-selectedTargetTime(isel));
    summaryRows(isel).NInitialAzimuths = numel(sdatAzimuth_deg);
    summaryRows(isel).Nvalid440_beforeCommonMask = nValidPerWL(1);
    summaryRows(isel).Nvalid675_beforeCommonMask = nValidPerWL(2);
    summaryRows(isel).Nvalid870_beforeCommonMask = nValidPerWL(3);
    summaryRows(isel).Nvalid1020_beforeCommonMask = nValidPerWL(4);
    summaryRows(isel).AllWL_atLeast10_beforeCommonMask = all(nValidPerWL >= 10);
    summaryRows(isel).AllWL_haveAll4AngularRanges = all(rangeOKperWL);
    summaryRows(isel).N_SDAT_common = sum(commonKeep);
    summaryRows(isel).Nremoved_byCommonMask = sum(~commonKeep);
    summaryRows(isel).RemovedAzimuth_deg = join(string(removedAz'), ';');
    summaryRows(isel).RetainedAzimuth_deg = join(string(retainedAz'), ';');
    summaryRows(isel).RetainedRAA_deg = join(string(retainedRAA'), ';');
    summaryRows(isel).Ncommon_2to6 = n2to6;
    summaryRows(isel).Ncommon_6to30 = n6to30;
    summaryRows(isel).Ncommon_30to80 = n30to80;
    summaryRows(isel).Ncommon_gt80 = ngt80;

    summaryRows(isel).SZA440_deg = sza4(1);
    summaryRows(isel).SZA675_deg = sza4(2);
    summaryRows(isel).SZA870_deg = sza4(3);
    summaryRows(isel).SZA1020_deg = sza4(4);
    summaryRows(isel).VZA440_deg = 180-sza4(1);
    summaryRows(isel).VZA675_deg = 180-sza4(2);
    summaryRows(isel).VZA870_deg = 180-sza4(3);
    summaryRows(isel).VZA1020_deg = 180-sza4(4);
    summaryRows(isel).RowTime440_UTC = rowTime4(1);
    summaryRows(isel).RowTime675_UTC = rowTime4(2);
    summaryRows(isel).RowTime870_UTC = rowTime4(3);
    summaryRows(isel).RowTime1020_UTC = rowTime4(4);

    % One row for every initial azimuth.
    for jang = 1:numel(sdatAzimuth_deg)
        ia = ia+1;
        angleRows(ia).TargetIndex = selectedTargetIndex(isel);
        angleRows(ia).Label = selectedLabel(isel);
        angleRows(ia).MatchedScanCenter_UTC = s.CenterTime;
        angleRows(ia).Azimuth_deg = sdatAzimuth_deg(jang);
        angleRows(ia).RAA_deg = sdatAzimuth_deg(jang)+180;
        angleRows(ia).RetainedInSDAT = commonKeep(jang);
        angleRows(ia).Valid440 = valid4(1,jang);
        angleRows(ia).Valid675 = valid4(2,jang);
        angleRows(ia).Valid870 = valid4(3,jang);
        angleRows(ia).Valid1020 = valid4(4,jang);
        angleRows(ia).RawRadiance440 = rad4(1,jang);
        angleRows(ia).RawRadiance675 = rad4(2,jang);
        angleRows(ia).RawRadiance870 = rad4(3,jang);
        angleRows(ia).RawRadiance1020 = rad4(4,jang);
    end
end

summaryTable = struct2table(summaryRows);
angleTable = struct2table(angleRows);

%% --------------------------- DISPLAY ---------------------------------

fprintf('\n================ SDAT CROSS-CHECK SUMMARY ================\n');
disp(summaryTable(:, { ...
    'Label','RequestedTime_UTC','MatchedScanCenter_UTC','NInitialAzimuths', ...
    'Nvalid440_beforeCommonMask','Nvalid675_beforeCommonMask', ...
    'Nvalid870_beforeCommonMask','Nvalid1020_beforeCommonMask', ...
    'N_SDAT_common','RemovedAzimuth_deg', ...
    'SZA440_deg','SZA675_deg','SZA870_deg','SZA1020_deg'}));

fprintf('\nKey fields for comparison with an independently generated SDAT:\n');
fprintf('  N_SDAT_common          : number of angular measurements retained.\n');
fprintf('  RetainedAzimuth_deg    : azimuths after the common 4-WL mask.\n');
fprintf('  RetainedRAA_deg        : azimuth + 180.\n');
fprintf('  SZAxxx_deg / VZAxxx_deg: wavelength-specific geometry.\n');
fprintf('  angle-level CSV        : identifies the wavelength causing removal.\n');
fprintf('  Scattering angle is not used in this SDAT cross-check.\n');

%% ---------------------------- SAVE -----------------------------------

summaryCsv = fullfile(outDir,'almucantar_sdat_crosscheck_summary.csv');
anglesCsv = fullfile(outDir,'almucantar_sdat_crosscheck_angles.csv');
writetable(summaryTable,summaryCsv);
writetable(angleTable,anglesCsv);

fprintf('\nSaved:\n  %s\n  %s\n',summaryCsv,anglesCsv);

%% ========================== LOCAL FUNCTION ============================

function [almFile, cleanupObj] = getAlmFile(inputFile)
% Accept an extracted .alm file or a ZIP containing one .alm file.
    inputFile = char(inputFile);
    cleanupObj = [];

    if endsWith(lower(inputFile),'.zip')
        tmpDir = tempname;
        mkdir(tmpDir);
        unzip(inputFile,tmpDir);
        d = dir(fullfile(tmpDir,'**','*.alm'));
        if isempty(d)
            rmdir(tmpDir,'s');
            error('ZIP does not contain an .alm file: %s',inputFile);
        end
        if numel(d) > 1
            rmdir(tmpDir,'s');
            error('ZIP contains more than one .alm file; extract and select the intended file.');
        end
        almFile = fullfile(d(1).folder,d(1).name);
        cleanupObj = onCleanup(@() rmdir(tmpDir,'s'));
    elseif endsWith(lower(inputFile),'.alm')
        almFile = inputFile;
    else
        error('inputFile must be an .alm file or a .zip containing one .alm file.');
    end
end
