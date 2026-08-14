%% grasp_avp_normalization_check.m
% Check normalization of the aerosol vertical profile (AVP) printed by
% GRASP classic inversion output files.
%
% The script reads one or more *_inversion_output.txt files and calculates:
%   1) integral of the 60 retrieved AVP points;
%   2) contribution of the extra ground point printed by GRASP;
%   3) integral of the complete printed AVP (including ground extension);
%   4) missing fraction required for the integral to equal 1;
%   5) the constant-continuation depth that would supply the missing area
%      if the uppermost printed AVP value were held constant above z_top;
%   6) the corresponding implied upper altitude;
%   7) AOD1064 * integral(printed AVP), for comparison with an AOD obtained
%      by integrating extinction reconstructed as AOD1064 * AVP.
%
% IMPORTANT INTERPRETATION:
% The inferred continuation depth / upper altitude is a DIAGNOSTIC only.
% Agreement between cases is evidence for a fixed upper-boundary treatment,
% but does not by itself prove how GRASP implements the profile internally.
% Confirmation from the GRASP developers/source code is still required.
%
% No special MATLAB toolboxes are required.

clear; clc;

%% ------------------------ USER SETTINGS -------------------------------

% Option A: automatically analyse all matching files in the current folder.
files = dir('*_inversion_output.txt');

% Option B: use an explicit list instead (uncomment and edit as needed).
% fileNames = {
%     '20240713T143509_inversion_output.txt'
%     '20240713T145353_inversion_output.txt'
%     '20240904T134749_inversion_output.txt'
%     '20240923T071917_inversion_output.txt'
%     };
% files = struct('name', fileNames);

outputCsv = 'grasp_avp_normalization_summary.csv';

if isempty(files)
    error('No *_inversion_output.txt files found in the current folder.');
end

%% ------------------------- CALCULATIONS -------------------------------

nFiles = numel(files);

File = strings(nFiles,1);
DateUTC = strings(nFiles,1);
TimeUTC = strings(nFiles,1);
AOD1064 = nan(nFiles,1);
AVPIntegralRetrieved = nan(nFiles,1);
GroundExtensionContribution = nan(nFiles,1);
AVPIntegralPrinted = nan(nFiles,1);
MissingToUnity = nan(nFiles,1);
TopAltitude_mASL = nan(nFiles,1);
TopAVP_mInv = nan(nFiles,1);
ImpliedContinuationDepth_m = nan(nFiles,1);
ImpliedUpperAltitude_mASL = nan(nFiles,1);
AOD_from_AODxPrintedAVP = nan(nFiles,1);
AODRatio_fromPrintedAVP = nan(nFiles,1);

for i = 1:nFiles
    fileName = files(i).name;
    File(i) = string(fileName);
    txt = fileread(fileName);

    %% Date and time
    tok = regexp(txt, 'Date:\s*(\d{4}-\d{2}-\d{2})', 'tokens', 'once');
    if ~isempty(tok)
        DateUTC(i) = string(tok{1});
    end

    tok = regexp(txt, 'Time:\s*(\d{2}:\d{2}:\d{2})', 'tokens', 'once');
    if ~isempty(tok)
        TimeUTC(i) = string(tok{1});
    end

    %% Extract AVP block
    startMarker = 'Aerosol vertical profile [1/m] for Particle component 1';
    endMarker   = 'Angstrom exp';

    iStart = strfind(txt, startMarker);
    if isempty(iStart)
        error('AVP section not found in %s', fileName);
    end
    iStart = iStart(1) + length(startMarker);

    txtAfterStart = txt(iStart:end);
    iEndRel = strfind(txtAfterStart, endMarker);
    if isempty(iEndRel)
        error('End of AVP section not found in %s', fileName);
    end

    avpBlock = txtAfterStart(1:iEndRel(1)-1);
    lines = regexp(avpBlock, '\r\n|\n|\r', 'split');

    idx = [];
    z = [];
    avp = [];
    isGround = [];

    for j = 1:numel(lines)
        % Captures both normal lines, e.g.
        %   60    242.50   0.23186E-03
        % and the starred ground line, e.g.
        % * 61     70.00   0.23186E-03
        t = regexp(lines{j}, ...
            '^\s*(\*?)\s*(\d+)\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s*$', ...
            'tokens', 'once');

        if ~isempty(t)
            idx(end+1,1) = str2double(t{2}); %#ok<SAGROW>
            z(end+1,1) = str2double(t{3}); %#ok<SAGROW>
            avp(end+1,1) = str2double(t{4}); %#ok<SAGROW>
            isGround(end+1,1) = ~isempty(t{1}); %#ok<SAGROW>
        end
    end

    if numel(z) < 2
        error('Could not parse AVP values from %s', fileName);
    end

    % Sort by physical altitude because GRASP prints the profile from top
    % to bottom, while trapz should be applied with increasing z.
    [zAll, orderAll] = sort(z, 'ascend');
    avpAll = avp(orderAll);
    groundAll = logical(isGround(orderAll));

    % Complete printed profile, including the starred ground point.
    Iprinted = trapz(zAll, avpAll);

    % Retrieved points only (exclude the starred ground point).
    zRet = zAll(~groundAll);
    avpRet = avpAll(~groundAll);

    if numel(zRet) < 2
        error('Fewer than two retrieved AVP points in %s', fileName);
    end

    Iretrieved = trapz(zRet, avpRet);
    Iground = Iprinted - Iretrieved;

    % Uppermost retrieved point.
    [zTop, iTop] = max(zRet);
    avpTop = avpRet(iTop);

    missing = 1 - Iprinted;

    % If the final AVP value were continued constantly above zTop, this is
    % the extra vertical distance required for the total integral to reach 1.
    if missing >= 0 && avpTop > 0
        continuationDepth = missing / avpTop;
        impliedUpperAltitude = zTop + continuationDepth;
    else
        continuationDepth = NaN;
        impliedUpperAltitude = NaN;
    end

    %% Extract GRASP AOD at 1.064 um
    % Isolate the first AOD_Total section so later output sections cannot
    % accidentally be matched.
    aodStart = strfind(txt, 'Wavelength (um), AOD_Total');
    if isempty(aodStart)
        error('AOD_Total section not found in %s', fileName);
    end
    aodText = txt(aodStart(1):end);

    aodEnd = strfind(aodText, 'Wavelength (um), SSA_Total');
    if ~isempty(aodEnd)
        aodText = aodText(1:aodEnd(1)-1);
    end

    tAod = regexp(aodText, ...
        '^\s*1\.06400\s+([0-9.+\-Ee]+)\s*$', ...
        'tokens', 'once', 'lineanchors');

    if isempty(tAod)
        error('1.064 um AOD not found in %s', fileName);
    end

    aod1064 = str2double(tAod{1});

    %% Save case results
    AOD1064(i) = aod1064;
    AVPIntegralRetrieved(i) = Iretrieved;
    GroundExtensionContribution(i) = Iground;
    AVPIntegralPrinted(i) = Iprinted;
    MissingToUnity(i) = missing;
    TopAltitude_mASL(i) = zTop;
    TopAVP_mInv(i) = avpTop;
    ImpliedContinuationDepth_m(i) = continuationDepth;
    ImpliedUpperAltitude_mASL(i) = impliedUpperAltitude;
    AOD_from_AODxPrintedAVP(i) = aod1064 * Iprinted;
    AODRatio_fromPrintedAVP(i) = Iprinted;
end

%% ----------------------------- OUTPUT ---------------------------------

T = table(File, DateUTC, TimeUTC, AOD1064, ...
    AVPIntegralRetrieved, GroundExtensionContribution, AVPIntegralPrinted, ...
    MissingToUnity, TopAltitude_mASL, TopAVP_mInv, ...
    ImpliedContinuationDepth_m, ImpliedUpperAltitude_mASL, ...
    AOD_from_AODxPrintedAVP, AODRatio_fromPrintedAVP);

disp(T);
writetable(T, outputCsv);

validUpper = ImpliedUpperAltitude_mASL(isfinite(ImpliedUpperAltitude_mASL));

fprintf('\n');
fprintf('Saved summary: %s\n', outputCsv);

if ~isempty(validUpper)
    fprintf('Implied constant-continuation upper altitude across cases:\n');
    fprintf('  mean   = %.3f m a.s.l.\n', mean(validUpper));
    fprintf('  median = %.3f m a.s.l.\n', median(validUpper));
    fprintf('  min    = %.3f m a.s.l.\n', min(validUpper));
    fprintf('  max    = %.3f m a.s.l.\n', max(validUpper));
    fprintf('  range  = %.3f m\n', max(validUpper)-min(validUpper));
end

%% ----------------------- EXPECTED CROSS-CHECK -------------------------
% For the four test files discussed in August 2026, values should be close
% to the following (small differences can arise only from text rounding):
%
% 2024-07-13 14:35:09
%   printed AVP integral        ~ 0.961187
%   top AVP                    ~ 2.3426e-6 1/m
%   implied upper altitude     ~ 23432.4 m a.s.l.
%   AOD1064 * AVP integral     ~ 0.066138
%
% 2024-07-13 14:53:53
%   printed AVP integral        ~ 0.973122
%   top AVP                    ~ 1.6224e-6 1/m
%   implied upper altitude     ~ 23430.9 m a.s.l.
%   AOD1064 * AVP integral     ~ 0.069070
%
% 2024-09-04 13:47:49
%   printed AVP integral        ~ 0.674835
%   top AVP                    ~ 1.9626e-5 1/m
%   implied upper altitude     ~ 23432.0 m a.s.l.
%   AOD1064 * AVP integral     ~ 0.161826
%
% 2024-09-23 07:19:17
%   printed AVP integral        ~ 0.717833
%   top AVP                    ~ 1.7031e-5 1/m
%   implied upper altitude     ~ 23431.7 m a.s.l.
%   AOD1064 * AVP integral     ~ 0.041922
