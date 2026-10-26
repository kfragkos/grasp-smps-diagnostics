function [summary, angles] = almucantar_pair_comparison(inputFile)
% ALMUCANTAR_PAIR_COMPARISON  Reproducible diagnostics of paired scans.
%
% Reads a full AERONET Version 3 RAW almucantar (.alm) file (or a ZIP
% containing one .alm file), matches FOUR exact requested 2024 scan centres
% and reports raw/retained geometry and branch asymmetry.
%
% Screening: match the 30 outward azimuths (+2 to +180 deg, retaining
% just one +6) with their return-branch counterparts at -azimuth;
% reject a position if either radiance is missing/nonpositive or if
% 200*abs(Iplus-Iminus)/(Iplus+Iminus)>20% at ANY of the four wavelengths.
% The common 4-wavelength mask is applied at all wavelengths. Branch
% symmetry does not depend on radiance normalization (one factor per
% wavelength multiplies both branches). Scattering angles are used ONLY
% for diagnostic plots/coverage, not for constructing the SDAT mask.
%
% Outputs:
%   results/almucantar_case_pair_summary.csv
%   results/almucantar_case_pair_angles.csv
%   results/almucantar_case_pair_YYYYMMDD.png (when scans are present)
%
% This is an independent implementation of the reported photometer
% screening, not a byte-for-byte SDAT generator. Compare retained RAA
% against the actual SDAT before treating the masks as identical.
%
% Usage:
%   almucantar_pair_comparison('20240101_20241231_Magurele_Inoe.alm')
%   almucantar_pair_comparison('annual_almucantar.zip')
%
% Also works with shorter .alm files: absent requested scans are marked
% MISSING, not replaced by unrelated neighbouring scans.
% MATLAB R2020a+; no extra toolboxes.

if nargin < 1 || isempty(inputFile)
    inputFile = '20240101_20241231_Magurele_Inoe.alm';
end
wl = [440 675 870 1020];
targets = [
    datetime(2024,5,26,14,26,03,'TimeZone','UTC')
    datetime(2024,5,26,16,00,28,'TimeZone','UTC')
    datetime(2024,9,23,7,19,17,'TimeZone','UTC')
    datetime(2024,9,23,13,40,58,'TimeZone','UTC')
];
labels = ["26May_1426";"26May_1600";"23Sep_0719";"23Sep_1340"];
tolMinutes = 5;
branchThreshold = 20;
outputDir = 'results';
if ~isfolder(outputDir), mkdir(outputDir); end

[filename, cleanTmp] = getAlmFile(inputFile); %#ok<NASGU>
lines = readlines(filename);
hi = find(startsWith(lines,"AERONET_Site,Date("),1);
if isempty(hi), error('AERONET .alm header not found.'); end
h = split(lines(hi),',');
if h(end) == "", h(end) = []; end
nMeta = 11;
n = (numel(h)-nMeta)/2;
if n ~= floor(n), error('Unexpected radiance/scattering columns.'); end
azAll = str2double(h(nMeta+1:nMeta+n));
azAll = azAll(:);
idx180 = find(azAll==180,1);
positive = [];
negative = [];
for k = 1:idx180
    a = azAll(k);
    if a < 2 || any(abs(azAll(positive)-a)<1e-8), continue; end
    neg = find(azAll == -a & (1:n)' > idx180,1,'first');
    if ~isempty(neg)
        positive(end+1,1) = k; %#ok<AGROW>
        negative(end+1,1) = neg; %#ok<AGROW>
    end
end
az = azAll(positive);
if numel(az) ~= 30
    error('Expected 30 paired azimuth positions; found %d.', numel(az));
end

r = struct('Time',{},'WL',{},'SZA',{},'Rad',{},'Scat',{});
for k=hi+1:numel(lines)
    f=split(strtrim(lines(k)),',');
    if numel(f)<nMeta+2*n, continue; end
    thisWL = str2double(f(5));
    if ~ismember(thisWL,wl), continue; end
    try
        tm = datetime(f(2)+" "+f(3), ...
            'InputFormat','dd:MM:yyyy HH:mm:ss','TimeZone','UTC');
    catch
        continue;
    end
    a = numel(r)+1;
    r(a).Time=tm;
    r(a).WL=thisWL;
    r(a).SZA=str2double(f(7));
    r(a).Rad=str2double(f(nMeta+1:nMeta+n));
    r(a).Scat=str2double(f(nMeta+n+1:nMeta+2*n));
end
if isempty(r), error('No rows at the four GRASP wavelengths.'); end
[~,ord]=sort([r.Time]); r=r(ord);
scans = {};
s=1;
for k=2:numel(r)
    if minutes(r(k).Time-r(k-1).Time)>5
        scans{end+1}=s:k-1; %#ok<AGROW>
        s=k;
    end
end
scans{end+1}=s:numel(r);

% All target scans are represented, including those absent from input.
out=struct([]);
detail=struct([]);
d=0;
plotRecords=cell(4,1);
for it=1:4
    hit=0;
    best=Inf;
    matched=NaT(1,1,'TimeZone','UTC');
    chosen=[];
    for j=1:numel(scans)
        ii=scans{j};
        x=[r(ii).WL];
        if ~all(ismember(wl,x)), continue; end
        center=datetime(mean(posixtime([r(ii).Time])), ...
            'ConvertFrom','posixtime','TimeZone','UTC');
        err=abs(minutes(center-targets(it)));
        if err<best, best=err; hit=j; matched=center; end
    end
    found=(hit>0 && best<=tolMinutes);
    o.Label=labels(it); o.TargetTime_UTC=targets(it);
    o.MatchedScan_UTC=matched; o.Status="MISSING";
    o.Nretained=NaN; o.RemovedRAA_deg="";
    o.MeanSZA_deg=NaN;
    o.MaxRawScattering_deg=NaN; o.MaxRetainedScattering_deg=NaN;
    o.MeanRawBranchDiff_gt6_pct=NaN; o.MeanRetainedBranchDiff_gt6_pct=NaN;
    o.MaxRawBranchDiff_gt6_pct=NaN; o.MaxRetainedBranchDiff_gt6_pct=NaN;
    if ~found
        o.MatchedScan_UTC=NaT(1,1,'TimeZone','UTC');
        out=[out;o]; %#ok<AGROW>
        fprintf('%s: missing exact scan in input file.\n',labels(it));
        continue;
    end
    ii=scans{hit};
    radPlus=nan(4,30); radMinus=nan(4,30);
    scat=nan(4,30); sza=nan(4,1);
    for iw=1:4
        cand=ii([r(ii).WL]==wl(iw));
        [~,q]=min(abs(minutes([r(cand).Time]-matched)));
        row=r(cand(q));
        radPlus(iw,:)=row.Rad(positive)';
        radMinus(iw,:)=row.Rad(negative)';
        scat(iw,:)=mean([abs(row.Scat(positive)'); ...
                          abs(row.Scat(negative)')],1);
        sza(iw)=row.SZA;
    end
    denom=radPlus+radMinus;
    good=isfinite(radPlus)&isfinite(radMinus)& ...
         radPlus>0&radMinus>0&denom>0;
    delta=200*abs(radPlus-radMinus)./denom;
    delta(~good)=NaN;
    keepWL=good&(delta<=branchThreshold);
    keepAll=all(keepWL,1);
    plotRecords{it}=struct('Scat',scat,'Diff',delta,'KeepAll',keepAll);
    outer=(az'>6);
    evalRaw=outer & all(isfinite(delta),1);
    rawVals=delta(:,evalRaw);
    keepOuter=outer & keepAll;
    keptVals=delta(:,keepOuter);
    o.Status="MATCHED";
    o.Nretained=sum(keepAll);
    o.RemovedRAA_deg=join(string((az(~keepAll')+180)'),';');
    o.MeanSZA_deg=mean(sza);
    rawScat=scat(isfinite(scat) & scat>0);
    retainedScat=scat(:,keepAll);
    retainedScat=retainedScat(isfinite(retainedScat) & retainedScat>0);
    if ~isempty(rawScat), o.MaxRawScattering_deg=max(rawScat); end
    if ~isempty(retainedScat), o.MaxRetainedScattering_deg=max(retainedScat); end
    if ~isempty(rawVals)
        o.MeanRawBranchDiff_gt6_pct=mean(rawVals(:));
        o.MaxRawBranchDiff_gt6_pct=max(rawVals(:));
    end
    if ~isempty(keptVals)
        o.MeanRetainedBranchDiff_gt6_pct=mean(keptVals(:));
        o.MaxRetainedBranchDiff_gt6_pct=max(keptVals(:));
    end
    out=[out;o]; %#ok<AGROW>
    for ia=1:30
        d=d+1;
        detail(d).Label=labels(it);
        detail(d).Azimuth_deg=az(ia);
        detail(d).RAA_deg=az(ia)+180;
        detail(d).Retained=keepAll(ia);
        for iw=1:4
            detail(d).(['BranchDifference_' num2str(wl(iw)) '_pct'])=delta(iw,ia);
            detail(d).(['ValidAfterBranchScreen_' num2str(wl(iw))])=keepWL(iw,ia);
            detail(d).(['ScatteringAngle_' num2str(wl(iw)) '_deg'])=scat(iw,ia);
        end
    end
end
summary=struct2table(out);
if isempty(detail)
    angles=table();
else
    angles=struct2table(detail);
end
writetable(summary,fullfile(outputDir,'almucantar_case_pair_summary.csv'));
if ~isempty(angles)
    writetable(angles,fullfile(outputDir,'almucantar_case_pair_angles.csv'));
end
disp(summary);
for day=1:2
    ids=(day-1)*2+(1:2);
    if all(cellfun(@isempty,plotRecords(ids))), continue; end
    fig=figure('Color','w','Position',[50 50 1150 700]);
    for iw=1:4
        subplot(2,2,iw);hold on;grid on;box on;
        for k=ids
            x=plotRecords{k};
            if isempty(x),continue;end
            [xx,sortIdx]=sort(x.Scat(iw,:));
            yy=x.Diff(iw,sortIdx);
            plot(xx,yy,'-o','MarkerSize',3,'DisplayName',char(labels(k)));
            dropped=~x.KeepAll(sortIdx);
            scatter(xx(dropped),yy(dropped),28,'x', ...
                'HandleVisibility','off');
        end
        yline(20,'k--','20%','HandleVisibility','off');
        xlabel('Actual scattering angle [deg]');
        ylabel('Raw branch difference [%]');
        title(sprintf('%d nm',wl(iw)));legend('Location','best');
    end
    if day==1, tag='20240526';else,tag='20240923';end
    exportgraphics(fig,fullfile(outputDir, ...
        ['almucantar_case_pair_' tag '.png']),'Resolution',180);
    close(fig);
end
end

function [file,tmpClean]=getAlmFile(input)
tmpClean=[];
if endsWith(lower(string(input)),".alm")
    file=char(input);
elseif endsWith(lower(string(input)),".zip")
    tmp=tempname;mkdir(tmp); unzip(char(input),tmp);
    f=dir(fullfile(tmp,'**','*.alm'));
    if numel(f)~=1,error('ZIP must contain exactly one .alm file.');end
    file=fullfile(f(1).folder,f(1).name);
    tmpClean=onCleanup(@()rmdir(tmp,'s'));
else
    error('Input must be an .alm or a ZIP containing one .alm.');
end
end
