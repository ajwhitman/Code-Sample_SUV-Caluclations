maxNumCompThreads(3);
tStart = tic;
fprintf('%+81s',sprintf('Started at %s\n',datestr(datetime('now'))))
%% Setting up, reading in, formatting
addpath('/home/jagust/slbaker/Code/PET/') %/CalcR.m
arda = '/home/jagust/arda/lblid/';
av = 'AV1451*';
fprintf('%-55s%+25s\n','Putting together the data table... ',datestr(datetime('now')))
T = readtable('/home/jagust/whitman/jnj_things/AV1451_MultiScan_Data.csv', ...
    'ReadVariableNames',1, ...'ReadRowNames',1,...
    'Format','%s %{MM/dd/yyyy}D %f %f %{MM/dd/yyyy}D %f %f %{MM/dd/yyyy}D %f %f');
T.Date1 = datetime(T.Date1,'Format','yyyy-MM-dd');
T.Date2 = datetime(T.Date2,'Format','yyyy-MM-dd');
T.Date3 = datetime(T.Date3,'Format','yyyy-MM-dd');
T.Exist3 = ~isnat(T.Date3);
nsubs = height(T);

%       nsubs=3; % for practice, so it doesn't take forever all the time
%       fprintf('(Only doing %d subjects)', nsubs);
%       T=head(T,nsubs);

%% Getting main scan filepath
fprintf('%-55s%+25s\n','Gathering filepaths and such... ',datestr(datetime('now')))
tmp = arrayfun(@(r,d)(dir(fullfile(arda,char(r),[av char(d)]))),T.Subject,T.Date1);
T.Scan1 = arrayfun(@(f)(fullfile(f.folder,f.name)),tmp,'UniformOutput',false);

tmp = arrayfun(@(r,d)(dir(fullfile(arda,char(r),[av char(d)]))),T.Subject,T.Date2);
T.Scan2 = arrayfun(@(f)(fullfile(f.folder,f.name)),tmp,'UniformOutput',false);

tmp = arrayfun(@(r,d)(dir(fullfile(arda,char(r),[av char(d)]))),T.Subject,T.Date3,'UniformOutput',false);
tmp(~cellfun('isempty',tmp)) = cellfun(@(f)(fullfile(f.folder,f.name)),tmp(~cellfun('isempty',tmp)),'UniformOutput',false);
tmp(cellfun('isempty',tmp)) = {''};
T.Scan3 = tmp;

%% Getting paths of files we'll need to read

% suvr_cereg.nii % mean_size_cereg.csv % raparc+aseg.nii % rwrCerebellum-SUIT.nii % inf_cerebellar_mask.nii

files = [{'suvr_cereg.*.nii'} {'mean_size_cereg.csv$'} {'raparc\+aseg.nii'} ...
    {'rwrCerebellum-SUIT.nii'} {'inf_cerebellar_mask.nii'}];

for i=1:length(files) %Define varname, spm_select filepaths
    fil = files{i};
    varname = matlab.lang.makeValidName(strtok(fil,'.'));
    
    a = spm_select('FPList',T.Scan1,['^' fil]);
    T = addvars(T,a,'NewVariableNames',[varname '1'],'After','Scan1');
    
    a = spm_select('FPList',T.Scan2,['^' fil]);
    T = addvars(T,a,'NewVariableNames',[varname '2'],'After','Scan2');
    
    a = spm_select('FPList',T.Scan3,['^' fil]);
    tmp = cellstr(repmat(' ',nsubs,1));
    tmp(T.Exist3) = cellstr(a);
    T = addvars(T,tmp,'NewVariableNames',[varname '3'],'After','Scan3');
end

%% Reading in meaninfcereg from csv
fprintf('%-55s%+25s\n','Reading in meaninfcereg... ',datestr(datetime('now')))
%lastly you need to change the suvr_cereg.nii back to scanner units.
%to do this, read in the mean infcereg scanner units in mean_size_cereg.csv.
%As an example, for B10-282/AV1451 from 2014, the meaninfcereg = 3119.6085

tmp = arrayfun(@(f)(textscan(fopen(char(f)),'%f',1,'HeaderLines',3)),cellstr(T.mean_size_cereg1));
T = addvars(T,cell2mat(tmp),'NewVariableNames','meaninfcereg1','After','mean_size_cereg1');

tmp = arrayfun(@(f)(textscan(fopen(char(f)),'%f',1,'HeaderLines',3)),cellstr(T.mean_size_cereg2));
T = addvars(T,cell2mat(tmp),'NewVariableNames','meaninfcereg2','After','mean_size_cereg2');

tmp = arrayfun(@(f)(fopen(char(f))),cellstr(T.mean_size_cereg3));
tmp(tmp>0) = cell2mat(arrayfun(@(f)(textscan(f,'%f',1,'HeaderLines',3)),tmp(tmp>0)));
T = addvars(T,tmp,'NewVariableNames','meaninfcereg3','After','mean_size_cereg3');

fclose('all');

fprintf('%-55s%+25s\n','  Done gathering data, time to do some calculations!',datestr(datetime('now')))

%% Do and save the SUV calculations
S = table2struct(T);

ticSUV = tic;
fprintf('%-55s%+25s\n','Calculating SUVs... (this is the long part)',datestr(datetime('now')))

S = arrayfun(@calculateSUVs, S, 'UniformOutput',false); %Works because S is an array of structs
S = cell2mat(S); %For some reason, the previous function spits out a cell array of structs? 

tocSUV = toc(ticSUV);
fprintf('  Done calculating SUVs! It took %02.0fh %02.0fm %04.1fs \n',...
    fix(tocSUV/3600),mod(fix(tocSUV/60),60),mod(tocSUV,60))
%% Gather data to calculate slopes

% Gather data; dates for x-axis, then subregions / masks
Sdates = mat2cell(transpose(datenum([S.Date1;S.Date2;S.Date3])),ones(1,nsubs)); %Collect all dates
Sdates = cellfun(@(s)(s(~isnan(s))),Sdates,'UniformOutput',false);              %Get rid of NaNs/nonexsistant 3rd scan

N=34;
subsuvs = mat2cell(transpose([S.SUV_Subregions1;S.SUV_Subregions2;S.SUV_Subregions3]),ones(1,nsubs),[(N*2) N]); %Collect subregion suvs
subsuvs = cellfun(@(s)(s(repmat(~all(isnan(s)),size(s)))),subsuvs,'UniformOutput',false); %get rid of nonexistant 3rd scans
subsuvs = cellfun(@(a,b)([a b]),subsuvs(:,1),subsuvs(:,2),'UniformOutput',false); %concatenate columns (1st&2nd scan + 3rd scan)
subsuvs = cellfun(@(S)(reshape(S,[N,(length(S)/N)])),subsuvs,'UniformOutput',false); %Reshape into 2-3 cols, 34 rows

masksuvs = mat2cell(transpose([S.SUV_InfCereg1;S.SUV_InfCereg2;S.SUV_InfCereg3]),ones(1,nsubs)); %Collect all mask suvs
masksuvs = cellfun(@(s)(s(~isnan(s))),masksuvs,'UniformOutput',false); %get rid of NaNs/nonexistand 3rd scan

whitesuvs = mat2cell(transpose([S.SUV_White1;S.SUV_White2;S.SUV_White3]),ones(1,nsubs)); %Collect all white matter suvs
whitesuvs = cellfun(@(s)(s(~isnan(s))),whitesuvs,'UniformOutput',false); %get rid of NaNs/nonexistand 3rd scan

%% Calculate and save slopes

fprintf('%-55s%+25s\n','Calculating slopes... ',datestr(datetime('now')))

subslopes = cellfun(@calculateSlopes,Sdates,subsuvs,'UniformOutput',false);
[S.SubregionSlopes_byDay] = subslopes{:};
tmp = mat2cell([subslopes{:}]*365.25,N,ones(1,nsubs)); %made annual
[S.SubregionSlopes] = tmp{:};

maskslopes = cellfun(@calculateSlopes,Sdates,masksuvs);
tmp = num2cell(maskslopes);
[S.MaskSlopes_byDay] = tmp{:};
tmp = num2cell(maskslopes*365.25); %made annual
[S.MaskSlopes] = tmp{:};
% S = setfield(S,{nsubs},'MaskSlopes_byDay',maskslopes);
% S = setfield(S,{nsubs},'MaskSlopes',maskslopes*365.25); %made annual

whiteslopes = cellfun(@calculateSlopes,Sdates,whitesuvs);
tmp = num2cell(whiteslopes);
[S.WhiteSlopes_byDay] = tmp{:};
tmp = num2cell(whiteslopes*365.25); %made annual
[S.WhiteSlopes] = tmp{:};

%% Save variables, calculate percent change

fprintf('%-55s%+25s\n','Calculating percent change... ',datestr(datetime('now')))

subregion_slopes = [S.SubregionSlopes];
inf_cerebellar_mask_slopes = reshape([S.MaskSlopes],[1 nsubs]);
white_matter_slopes = reshape([S.WhiteSlopes],[1 nsubs]);

mean_subregion_slopes = nanmean(subregion_slopes,2);
mean_inf_cerebellar_mask_slopes = nanmean(inf_cerebellar_mask_slopes);
mean_white_matter_slopes = nanmean(white_matter_slopes);

perc_change_subregions = subregion_slopes./[S.SUV_Subregions1];
tmp = mat2cell(perc_change_subregions,N,ones(1,nsubs));
[S.PercChangeSubregions] = tmp{:};

perc_change_masks = inf_cerebellar_mask_slopes./[S.SUV_InfCereg1];
tmp = num2cell(perc_change_masks);
[S.PercChangeMasks] = tmp{:};

perc_change_white = white_matter_slopes./[S.SUV_White1];
tmp = num2cell(perc_change_white);
[S.PercChangeMasks] = tmp{:};

fprintf('%-55s%+25s\n','  Done calculating slopes and percent change.',datestr(datetime('now')))

%% Save output

fprintf('%-55s%+25s\n','Saving results... ',datestr(datetime('now')))

save(sprintf('SUV_Calculations_%s.mat',datestr(datetime('today'),'yyyy-mm-dd')),...
    'T','S','Sdates','nsubs','N',...
    'subsuvs','masksuvs','whitesuvs',...
    'subregion_slopes','inf_cerebellar_mask_slopes','white_matter_slopes',...
    'mean_subregion_slopes','mean_inf_cerebellar_mask_slopes','mean_white_matter_slopes',...
    'perc_change_subregions','perc_change_masks','perc_change_white');
fprintf('%-55s%+25s\n','Done!',datestr(datetime('now')))
%%
tStop = toc(tStart);
fprintf('Whole program took %02.0fh %02.0fm %04.1fs \n',...
    fix(tStop/3600),mod(fix(tStop/60),60),mod(tStop,60))
%% Plots prep
sr_sd = nanstd(subregion_slopes,0,2);
out_subject = mode(fix(find(abs(subregion_slopes) > max(sr_sd)*3)/N+1));
out_index = 1:nsubs ~= out_subject;
%subregion_slopes_noOutSubj = subregion_slopes(:,1:nsubs ~= out_subject);
%perc_change_subregions_noOutSubj = perc_change_subregions(:,1:nsubs ~= out_subject);

% out_points = find(abs(subregion_slopes) > max(sr_sd)*3);
% subregion_slopes_noOutPoints=subregion_slopes; subregion_slopes_noOutPoints(out_points) = nan;
% perc_change_subregions_noOutPoints=perc_change_subregions; perc_change_subregions_noOutPoints(out_points)=nan;

%% Plots
close all;
%A=ceil(max(max(abs(subregion_slopes)))*10)/10; B=round(max(sr_sd)*3,1,'significant');
%A=inf; B=inf; Y=inf; Y=45;
%axA = [-A A 0 Y]; axB = [-B B 0 Y];
%axAinf = [-A A 0 inf]; axBinf = [-B B 0 inf];
p=0; w=500; h=800; b=(1290-h); pint=w;

figure('Position', [p b w h]); 
suptitle('Slopes (all subjects)');
s1=subplot(3,1,1); histogram(subregion_slopes); 
ylabel('Frequency');xlabel('All subregion slopes');
%axis tight; maxlim = max(abs(get(s1,'XLim')));
%set(s1,'XLim',[-maxlim maxlim]); axis 'auto y';
s2=subplot(3,1,2); histogram(inf_cerebellar_mask_slopes);
ylabel('Frequency');xlabel('Inf Cere Mask slopes'); 
%axis tight; maxlim = max(abs(get(s1,'XLim')));
%set(s2,'XLim',[-maxlim maxlim]); axis 'auto y';
s3=subplot(3,1,3); histogram(white_matter_slopes);
ylabel('Frequency');xlabel('White Matter slopes'); 
%axis tight; maxlim = max(abs(get(s3,'XLim')));
%set(s3,'XLim',[-maxlim maxlim]); axis 'auto y';
linkaxes([s1,s2,s3],'x')
p=p+pint;

figure('Position', [p b w h]);
suptitle('Slopes (excluding subj 28)')
s1=subplot(3,1,1); histogram(subregion_slopes(:,out_index));
ylabel('Frequency');xlabel('All subregion slopes');
%axis tight; maxlim = max(abs(get(s1,'XLim')));
%set(s1,'XLim',[-maxlim maxlim]); axis 'auto y'
s2=subplot(3,1,2); histogram(inf_cerebellar_mask_slopes(:,out_index));
ylabel('Frequency');xlabel('Inf Cere Mask slopes');
%axis tight; maxlim = max(abs(get(s2,'XLim')));
%set(s2,'XLim',[-maxlim maxlim]); axis 'auto y';
s3=subplot(3,1,3); histogram(white_matter_slopes(:,out_index));
ylabel('Frequency');xlabel('White Matter slopes');
%axis tight; maxlim = max(abs(get(s3,'XLim')));
%set(s3,'XLim',[-maxlim maxlim]); axis 'auto y';
linkaxes([s1,s2,s3],'x')
p=p+pint;

figure('Position', [p b w h]);
suptitle('Percent change (all subjects)');
s1=subplot(3,1,1); histogram(perc_change_subregions);
ylabel('Frequency');xlabel('Percent change - All subregions');
s2=subplot(3,1,2); histogram(perc_change_masks);
ylabel('Frequency');xlabel('Percent change - Inf Cere Mask');
s3=subplot(3,1,3); histogram(perc_change_white);
ylabel('Frequency');xlabel('Percent change - White Matter');
linkaxes([s1,s2,s3],'x')
p=p+pint;

figure('Position', [p b w h]); 
suptitle('Percent change (excluding subj 28)');
s1=subplot(3,1,1); histogram(perc_change_subregions(:,out_index));
ylabel('Frequency');xlabel('Percent change - All subregions');
s2=subplot(3,1,2); histogram(perc_change_masks(:,out_index));
ylabel('Frequency');xlabel('Percent change - Inf Cere Mask');
s3=subplot(3,1,3); histogram(perc_change_white(:,out_index));
ylabel('Frequency');xlabel('Percent change - White Matter');
linkaxes([s1,s2,s3],'x')
p=p+pint;

figure;
plot(inf_cerebellar_mask_slopes(out_index),white_matter_slopes(out_index),'*')
axis([-.3 .3 -.3 .3]); axis square
hold on
plot([-.3 .3],[-.3 .3],'r-')
tmp=CalcR(inf_cerebellar_mask_slopes(out_index),white_matter_slopes(out_index));
x=[-.3 .3]; y=x.*tmp.p(1)-tmp.p(2);
plot(x,y,'-','LineWidth',2)
xlabel(['Annual Change SUV InfCereg'])
ylabel(['Annual Change SUV ErodedWhite'])

%%
subrHists1 = figure; set(subrHists1,'Position',[0 400 2000 800]);
suptitle('Subregion histograms (all subjects)');
sax = zeros(N,1);
for i=1:N
    sax(i)=subplot(2,N/2,i); 
    histogram(subregion_slopes(i,:));
    xlabel(sprintf('SubReg %d',i))
    %axis(axA)
end
linkaxes(sax,'xy');

subrHists2 = figure; set(subrHists2,'Position',[100 300 2000 800]);
suptitle('Subregion histograms (excluding subj 28)');
sax = zeros(N,1);
for i=1:N
    sax(i)=subplot(2,N/2,i); 
    histogram(subregion_slopes(i,out_index));
    xlabel(sprintf('SubReg %d',i))
    %axis(axB)
end
linkaxes(sax,'xy');

fprintf('Plots done\n')
%% Funcion to calculate slopes
function slopes = calculateSlopes(dates,suvarray) %(A) %(dates,suvarray)

%%% Make it Anual (/365.25) --edit: made annual after function call

%Pick only regions with multiple points 
goodones = sum(~isnan(suvarray),2)>1;

%turn array into single-column of ordered pairs
[N,~] = size(suvarray);
suvcell = mat2cell(suvarray,ones(1,N));

slopes = repmat(struct('r',NaN,'r2',NaN,'p',[NaN,NaN]),N,1); %default to NaN if no calculation
slopes(goodones) = cellfun(@(s)(CalcR(dates,s)),suvcell(goodones));
%Returns struct array with CalcR output; 
%Extract just the slope
slopes = cell2mat(transpose({slopes.p}));
slopes = slopes(:,1);
%slopes = slopes*365.25; %convert to years
end

%% Function to do all the work (calculating SUVs), that I can call from arrayfun b/c it's faster
function S = calculateSUVs(S)
fprintf('\tWorking on subject %s\n',S.Subject)
N = 34; %max(rwrC) = 34

    %for s=1:length(S) %height(T) %
    %Scan1:
        weight_dose = (S.Weight1*10^3)/(S.Dose1*37); % kg/mCi -> g/MBq
    
        %load rwrC, suvr_cereg, raparc+aseg; reshape to be 1D
        rwrC = spm_vol(S.rwrCerebellum_SUIT1);
        rwrC = spm_read_vols(rwrC);
        rwrC = reshape(rwrC,256*256*256,1);

        suvr = spm_vol(S.suvr_cereg1);
        suvr = spm_read_vols(suvr);
        suvr = reshape(suvr,256*256*256,1);

        raparc = spm_vol(S.raparc__aseg1);
        raparc = spm_read_vols(raparc); 
        raparc = reshape(raparc,256*256*256,1);
        
        
        suvarray1 = zeros(N,1);
        for n=1:N
            goodvox = rwrC==n & (raparc==8 | raparc==47); %find(rwrC==n & (raparc==8 | raparc==47));
            meansuvr = mean(suvr(goodvox));
            suv = (meansuvr*S.meaninfcereg1*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
            suvarray1(n) = suv;
        end
        
        %load inf_cerebellar_mask; reshape to be 1D
        infC = spm_vol(S.inf_cerebellar_mask1);
        infC = spm_read_vols(infC);
        infC = reshape(infC,256*256*256,1);
        goodvox = infC==1 & (raparc==8 | raparc==47);
        meansuvr = mean(suvr(goodvox));
        suv_infC1 = (meansuvr*S.meaninfcereg1*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
        
        %white matter suv
        rwhite=zeros(256.^3,1);
        ind=(raparc==2 | raparc==41);
        rwhite(ind)=1;
        white=reshape(rwhite,256,256,256);
        swhite=zeros(256,256,256);
        spm_smooth(white,swhite,[8 8 8]);
        rswhite=reshape(swhite,256.^3,1);
        indw=(rswhite>0.8 & rwhite>.5);
        meansuvr = mean(suvr(indw));
        suv_white1 = (meansuvr*S.meaninfcereg1*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
        

    %Scan2:
        weight_dose = (S.Weight2*10^3)/(S.Dose2*37); % kg/mCi -> g/MBq
        
        %load rwrC, suvr_cereg, raparc+aseg; reshape to be 1D
        rwrC = spm_vol(S.rwrCerebellum_SUIT2);
        rwrC = spm_read_vols(rwrC);
        rwrC = reshape(rwrC,256*256*256,1);

        suvr = spm_vol(S.suvr_cereg2);
        suvr = spm_read_vols(suvr);
        suvr = reshape(suvr,256*256*256,1);

        raparc = spm_vol(S.raparc__aseg2);
        raparc = spm_read_vols(raparc);
        raparc = reshape(raparc,256*256*256,1);

        suvarray2 = zeros(N,1);
        for n=1:N
            goodvox = rwrC==n & (raparc==8 | raparc==47); %find(rwrC==n & (raparc==8 | raparc==47));
            meansuvr = mean(suvr(goodvox));
            suv = (meansuvr*S.meaninfcereg2*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
            suvarray2(n) = suv;
        end
        
        %load inf_cerebellar_mask; reshape to be 1D
        infC = spm_vol(S.inf_cerebellar_mask2);
        infC = spm_read_vols(infC);
        infC = reshape(infC,256*256*256,1);
        goodvox = infC==1 & (raparc==8 | raparc==47);
        meansuvr = mean(suvr(goodvox));
        suv_infC2 = (meansuvr*S.meaninfcereg2*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
        
        %white matter suv
        rwhite=zeros(256.^3,1);
        ind=(raparc==2 | raparc==41);
        rwhite(ind)=1;
        white=reshape(rwhite,256,256,256);
        swhite=zeros(256,256,256);
        spm_smooth(white,swhite,[8 8 8]);
        rswhite=reshape(swhite,256.^3,1);
        indw=(rswhite>0.8 & rwhite>.5);
        meansuvr = mean(suvr(indw));
        suv_white2 = (meansuvr*S.meaninfcereg2*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
    
    
    %Scan3:
    suvarray3 = NaN(N,1);
    suv_infC3 = NaN;
    suv_white3 = NaN;
    if S.Exist3
            weight_dose = (S.Weight3*10^3)/(S.Dose3*37); % kg/mCi -> g/MBq
            
            %load rwrC, suvr_cereg, raparc+aseg; reshape to be 1D
            rwrC = spm_vol(S.rwrCerebellum_SUIT3);
            rwrC = spm_read_vols(rwrC);
            rwrC = reshape(rwrC,256*256*256,1);
            
            suvr = spm_vol(S.suvr_cereg3);
            suvr = spm_read_vols(suvr);
            suvr = reshape(suvr,256*256*256,1);
            
            raparc = spm_vol(S.raparc__aseg3);
            raparc = spm_read_vols(raparc);
            raparc = reshape(raparc,256*256*256,1);
            
            for n=1:N
                goodvox = rwrC==n & (raparc==8 | raparc==47); %find(rwrC==n & (raparc==8 | raparc==47));
                meansuvr = mean(suvr(goodvox));
                suv = (meansuvr*S.meaninfcereg3*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
                suvarray3(n) = suv;
            end
            
            %load inf_cerebellar_mask; reshape to be 1D
            infC = spm_vol(S.inf_cerebellar_mask3);
            infC = spm_read_vols(infC);
            infC = reshape(infC,256*256*256,1);
            goodvox = infC==1 & (raparc==8 | raparc==47);
            meansuvr = mean(suvr(goodvox));
            suv_infC3 = (meansuvr*S.meaninfcereg3*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
            
            %white matter suv
            rwhite=zeros(256.^3,1);
            ind=(raparc==2 | raparc==41);
            rwhite(ind)=1;
            white=reshape(rwhite,256,256,256);
            swhite=zeros(256,256,256);
            spm_smooth(white,swhite,[8 8 8]);
            rswhite=reshape(swhite,256.^3,1);
            indw=(rswhite>0.8 & rwhite>.5);
            meansuvr = mean(suvr(indw));
            suv_white3 = (meansuvr*S.meaninfcereg3*10^-6)*weight_dose; %(Bq/mL -> MBq/mL)*(g/MBq)
    end
    
    S.SUV_Subregions1 = suvarray1;
    S.SUV_InfCereg1   = suv_infC1;
    S.SUV_White1      = suv_white1;
    
    S.SUV_Subregions2 = suvarray2;
    S.SUV_InfCereg2   = suv_infC2;
    S.SUV_White2      = suv_white2;
    
    S.SUV_Subregions3 = suvarray3;
    S.SUV_InfCereg3   = suv_infC3;
    S.SUV_White3      = suv_white3;
    
end
