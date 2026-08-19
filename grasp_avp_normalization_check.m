%% grasp_avp_normalization_check.m
% Check normalization of the aerosol vertical profile (AVP) printed by
% GRASP classic inversion output files.
%
% For a LUT aerosol profile, GRASP adds boundary points at the ground and
% at the model top. Below the lowest LUT altitude the lowest value is held
% constant. Above the highest LUT altitude the profile is connected
% linearly to zero at the model top. For the configuration considered here,
% the top of the atmosphere is 40 km.
%
% The script reads one or more *_inversion_output.txt files and calculates:
%   1) integral of the retrieved AVP points;
%   2) contribution of the additional ground point printed by GRASP;
%   3) integral of the complete printed AVP up to the highest retrieved
%      altitude;
%   4) missing fraction required for the full normalized profile to equal 1;
%   5) expected upper-triangle contribution for a linear decrease to zero
%      at the configured model top;
%   6) model-top altitude inferred independently from the printed AVP and
%      the missing normalized fraction;
%   7) normalization closure after adding the upper triangle;
%   8) corresponding 1064-nm AOD contributions.
%
% No special MATLAB toolboxes are required.

clear; clc;

%% ------------------------ USER SETTINGS -------------------------------

% Analyse all matching GRASP classic inversion outputs in the current folder.
files = dir('*_inversion_output.txt');

% Alternatively, provide an explicit list:
% fileNames = {
%     '20240713T143509_inversion_output.txt'
%     '20240713T145353_inversion_output.txt'
%     '20240904T134749_inversion_output.txt'
%     '20240923T071917_inversion_output.txt'
%     };
% files = struct('name', fileNames);

% Top of the GRASP model atmosphere for this configuration.
HMAX_mASL = 40000;

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
UpperTriangle40km = nan(nFiles,1);
InferredHMAX_mASL = nan(nFiles,1);
NormalizationClosure40km = nan(nFiles,1);
ClosureError40km = nan(nFiles,1);
AOD_PrintRange = nan(nFiles,1);
AOD_UpperTriangle40km = nan(nFiles,1);
AOD_FullProfile40km = nan(nFiles,1);

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

    z = [];
    avp = [];
    isGround = [];

    for j = 1:numel(lines)
        % Captures both normal profile rows and the starred ground row.
        t = regexp(lines{j}, ...
            '^\s*(\*?)\s*(\d+)\s+([0-9.+\-Ee]+)\s+([0-9.+\-Ee]+)\s*$', ...
            'tokens', 'once');

        if ~isempty(t)
            z(end+1,1) = str2double(t{3}); %#ok<SAGROW>
            avp(end+1,1) = str2double(t{4}); %#ok<SAGROW>
            isGround(end+1,1) = ~isempty(t{1}); %#ok<SAGROW>
        end
    end

    if numel(z) < 2
        error('Could not parse AVP values from %s', fileName);
    end

    % GRASP prints top-to-bottom; sort by increasing physical altitude for trapz.
    [zAll, orderAll] = sort(z, 'ascend');
    avpAll = avp(orderAll);
    groundAll = logical(isGround(orderAll));

    % Complete printed profile, including the starred ground point.
    Iprinted = trapz(zAll, avpAll);

    % Retrieved points only.
    zRet = zAll(~groundAll);
    avpRet = avpAll(~groundAll);

    if numel(zRet) < 2
        error('Fewer than two retrieved AVP points in %s', fileName);
    end

    Iretrieved = trapz(zRet, avpRet);
    Iground = Iprinted - Iretrieved;

    [zTop, iTop] = max(zRet);
    avpTop = avpRet(iTop);
    missing = 1 - Iprinted;

    % GRASP LUT treatment above the highest retrieved altitude:
    % linear connection from AVP_top to zero at HMAX.
    if HMAX_mASL > zTop && avpTop >= 0
        upperTriangle = 0.5 * avpTop * (HMAX_mASL - zTop);
    else
        upperTriangle = NaN;
    end

    % Infer HMAX from the missing normalized area without assuming 40 km:
    % missing = 0.5 * AVP_top * (HMAX - zTop).
    if missing >= 0 && avpTop > 0
        inferredHmax = zTop + 2 * missing / avpTop;
    else
        inferredHmax = NaN;
    end

    closure = Iprinted + upperTriangle;
    closureError = closure - 1;

    %% Extract GRASP AOD at 1.064 um
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

    %% Save results
    AOD1064(i) = aod1064;
    AVPIntegralRetrieved(i) = Iretrieved;
    GroundExtensionContribution(i) = Iground;
    AVPIntegralPrinted(i) = Iprinted;
    MissingToUnity(i) = missing;
    TopAltitude_mASL(i) = zTop;
    TopAVP_mInv(i) = avpTop;
    UpperTriangle40km(i) = upperTriangle;
    InferredHMAX_mASL(i) = inferredHmax;
    NormalizationClosure40km(i) = closure;
    ClosureError40km(i) = closureError;
    AOD_PrintRange(i) = aod1064 * Iprinted;
    AOD_UpperTriangle40km(i) = aod1064 * upperTriangle;
    AOD_FullProfile40km(i) = aod1064 * closure;
end

%% ----------------------------- OUTPUT ---------------------------------

T = table(File, DateUTC, TimeUTC, AOD1064, ...
    AVPIntegralRetrieved, GroundExtensionContribution, AVPIntegralPrinted, ...
    MissingToUnity, TopAltitude_mASL, TopAVP_mInv, UpperTriangle40km, ...
    InferredHMAX_mASL, NormalizationClosure40km, ClosureError40km, ...
    AOD_PrintRange, AOD_UpperTriangle40km, AOD_FullProfile40km);

disp(T);
writetable(T, outputCsv);

validH = InferredHMAX_mASL(isfinite(InferredHMAX_mASL));

fprintf('\nSaved summary: %s\n', outputCsv);
fprintf('Configured model top: %.3f m a.s.l.\n', HMAX_mASL);

if ~isempty(validH)
    fprintf('Inferred model-top altitude across cases:\n');
    fprintf('  mean   = %.3f m a.s.l.\n', mean(validH));
    fprintf('  median = %.3f m a.s.l.\n', median(validH));
    fprintf('  min    = %.3f m a.s.l.\n', min(validH));
    fprintf('  max    = %.3f m a.s.l.\n', max(validH));
    fprintf('  range  = %.3f m\n', max(validH)-min(validH));
end

fprintf('\nInterpretation:\n');
fprintf('  AVPIntegralPrinted integrates only the printed profile range.\n');
fprintf('  UpperTriangle40km adds the linear extension from AVP_top to zero at 40 km.\n');
fprintf('  NormalizationClosure40km should be close to 1.\n');
fprintf('  InferredHMAX_mASL should be close to 40000 m if the printed profile is consistent with this treatment.\n');
