%% Find Tpmp, Tps, and Tpps, solve for H-k
%  Modified 7/18/2019
% This code performs H-k and H-k-Vp stacking and was developed to be
%    interactive. It should show both waveform and parallel plots to allow
%    the user to define: 1) whether gridsearch parameters need to be
%    broadened (Important to obtain useful errors), 2) whether cluster
%    analysis needs to be run to obtain appropriate solution.
%
% FUNCTION NECESSARY TO RUN
%   HkVp_times.m: Predicts arrivals times for H-k and H-k-Vp stacking
%   plot_HkVp.m: Makes plots of waveforms and H-k(-Vp) solutions
close all; clear all

% Folder where result matrices reside
datafolder = sprintf('../RESULTS_INTER');

% filename containing result matrices to make figures from
fid2 = fopen('results.txt','r');
Input = textscan(fid2,'%s %f %f');
results = Input{1};
STALOCS = [Input{2} Input{3}];
fclose(fid2);

% Toggles
w = [0.4 0.2 0.4] % weighting in order for Ps, 3P1S, and Pmp
wtog = 0; % normalize ACC weighting by number of events in bin...
          % 0 if unweighted, 1 if weighted by # traces in bin
ptog = 0; % If ptog = 1, amplitudes of opposite sign predicted along moveout 
          % curve are not summed for stack. If 0, all amplitudes along moveout
          % are stacked (can lead to destructive amplitude summing)
plotme = 0; % Whether to make and output plots (takes LONG time)

% Do stacking on Linear or PW stack
toggle = 0; % 0 for linear, 1 for pw
zonly = 1; % 1: only use Z autocorrelation amplitude data, 2: use Z and R

% SEARCH PARAMETERS
step_H = 0.2;
H = 25:step_H:60;

step_Vp = 0.02;
Vp = 5.6:step_Vp:7.2;

step_VpVs = 0.005;
VpVs = 1.65:step_VpVs:1.95;

SP = [min(H) max(H) step_H min(Vp) max(Vp) step_Vp min(VpVs) max(VpVs) step_VpVs];

% Assumed Vp for H-k stack comparison
AssVp = 6.2;

% RESULTS MATFILE
str1 = sprintf('ALL_RESULTS_inter_%0.2f_%0.2f_%0.2f',w(1),w(2),w(3));
if exist(strcat(str1,'.mat')) == 2
    load(strcat(str1,'.mat'))
else
    ALL_RESULT = zeros(length(results),25);
    SP = bsxfun(@times,ones(length(results),9),SP);
    clustvec = ones(length(results),1);
    ALL_RESULT = [ALL_RESULT SP clustvec];
    clear SP clustvec
end

%% DONT NEED TO CHANGE BELOW HERE
% for MATLAB's stupid plotting purposes only
%[VP,VPVS,HH1] = meshgrid(Vp,VpVs,H);
%[VPVS,HH2] = meshgrid(VpVs,H);
%H_VP_VPVS_mat = zeros(size(VP));
%H_VPVS_mat = zeros(size(VPVS));

%[XX,Z] = meshgrid(H,Vs);
if wtog == 1
    disp('Performing weighted H-k-Vp analysis');
else
    disp('Performing H-k-Vp analysis');
end

if ptog == 1
    disp('Eliminating amplitudes of opposite sign from stacks');
else
    disp('All amplitude info along moveout retained');
end


%% Make table of Tpmp, Tps, and Tpps travel times
% This makes code less robust as rayp below needs to match rayp
% centerpoints for bins, and must be done manually
rayp = [0.0402, 0.0446, 0.0490, 0.0534, 0.0578, 0.0622, 0.0666, 0.0710, 0.0754, 0.0798];

% Make vector of all station names to be analyzed
for i = 1:length(results)
    result = results(i);
    tmp2 = strsplit(char(result),'_');
    stnm = tmp2(1,1);
    STALST(i) = stnm;
end
STALST = STALST';

%% Run Stacking
for i = 1%:length(results) %%[ 3 5 6 7 12 25 38 52 54 ] %
    disp(sprintf('Loading Results for %s',results{i}));
    result = results(i);
    infile = strcat('./',datafolder,'/',result);
    load(sprintf('%s',infile{1}));
    stnm = STALST(i);

    %tmp2 = strsplit(char(result),'_');
    %stnm = tmp2(1,1)
    %STALST(i) = stnm;
    % Plotting and stacking info (0 to use linear stack amplitudes, 1 to 
    %use PW stack amplitudes
    if toggle == 0
        COL = 5;
        fname = strcat('./',datafolder,'/',char(stnm),'_HVpk_LIN.pdf');
        scale = 0.005;
    elseif toggle == 1
        COL = 6;
        fname = strcat('./',datafolder,'/',char(stnm),'_HVpk_PW.pdf');
        scale = 0.008;
    else
        disp('Stacking not specified. Using linear stacking')
        COL = 5;
        fname = strcat('./',datafolder,'/',char(stnm),'_HVpk_LIN.pdf');
        scale = 0.005;
    end
    % Name for output of gridsearch result file
    fname2 = strcat('./',datafolder,'/',char(stnm),'_terninfo.txt');
    fname21 = strcat('./',datafolder,'/',char(stnm),'_terninfo.mat');
    fname3 = strcat('./',datafolder,'/',char(stnm),'_BHinfo.txt');
    fname4 = strcat('./',datafolder,'/',char(stnm),'_clustsolns.txt');


    GOOD = 0;
    % Until good solution is declared
    while GOOD == 0
        tic
        SP = ALL_RESULT(i,26:34);
        w = ALL_RESULT(i,36:38);
        [tern_info, HkVp_RESULT, HkVp_good, Hk_good, std_err_HkVp, std_err_Hk] = HkVp(stack_ACCZ,stack_ACCR,stack_RFs,SP,rayp,zonly,AssVp,wtog,ptog,w,COL,fname2,fname21,fname3);

        ALL_RESULT(i,1)=HkVp_RESULT(1);ALL_RESULT(i,2)=HkVp_RESULT(2);ALL_RESULT(i,3)=HkVp_RESULT(3);ALL_RESULT(i,4)=HkVp_RESULT(4);ALL_RESULT(i,5)=HkVp_RESULT(5);ALL_RESULT(i,6)=HkVp_RESULT(6);
        ALL_RESULT(i,7)=HkVp_RESULT(7);ALL_RESULT(i,8)=HkVp_RESULT(8);ALL_RESULT(i,9)=HkVp_RESULT(9);ALL_RESULT(i,10)=HkVp_RESULT(10);ALL_RESULT(i,11)=HkVp_RESULT(11);ALL_RESULT(i,12)=HkVp_RESULT(12);
        ALL_RESULT(i,13)=HkVp_RESULT(13);ALL_RESULT(i,14)=HkVp_RESULT(14);ALL_RESULT(i,15)=HkVp_RESULT(15);ALL_RESULT(i,16)=HkVp_RESULT(16);ALL_RESULT(i,17)=HkVp_RESULT(17);ALL_RESULT(i,18)=HkVp_RESULT(18);
        ALL_RESULT(i,19)=HkVp_RESULT(19);ALL_RESULT(i,20)=HkVp_RESULT(20);ALL_RESULT(i,21)=HkVp_RESULT(21);ALL_RESULT(i,22)=HkVp_RESULT(22);ALL_RESULT(i,23)=HkVp_RESULT(23);ALL_RESULT(i,24)=HkVp_RESULT(24);    
        ALL_RESULT(i,25)=HkVp_RESULT(25);ALL_RESULT(i,26)=HkVp_RESULT(26);ALL_RESULT(i,27)=HkVp_RESULT(27);ALL_RESULT(i,28)=HkVp_RESULT(28);ALL_RESULT(i,29)=HkVp_RESULT(29);ALL_RESULT(i,30)=HkVp_RESULT(30);
        ALL_RESULT(i,31)=HkVp_RESULT(31);ALL_RESULT(i,32)=HkVp_RESULT(32);ALL_RESULT(i,33)=HkVp_RESULT(33);ALL_RESULT(i,34)=HkVp_RESULT(34);ALL_RESULT(i,36)=HkVp_RESULT(36);ALL_RESULT(i,37)=HkVp_RESULT(37);
        ALL_RESULT(i,38)=HkVp_RESULT(38);
%         if ALL_RESULT(i,35) == 0
%             n_clust = 1;
%             ALL_RESULT(i,35) = n_clust;
%         else
%             n_clust = ALL_RESULT(i,35);
%         end

        toc

        % PLOT  RESULTS AND CHECK SOLUTIONS
        % This function allows you to change parameters for HkVp/Hk
        % stacking if you want, but doesn't run HkVp stacking itself.
        % Cluster analysis, however, is run in this script.
        [ GOOD newSOLN] =plot_HkVp(fname,fname2,fname4,stack_ACCZ,stack_ACCR,stack_RFs,ALL_RESULT(i,:),std_err_HkVp,std_err_Hk,HkVp_good,Hk_good,COL,scale,plotme);
        ALL_RESULT(i,:) = newSOLN;

    end

    save(strcat(str1,'.mat'),'STALST','STALOCS','ALL_RESULT');

    %newSP = cell2mat(newSP);
    % WRITE IF (OR WHILE?) FOR NEW PARAMETER SEARCH
%     if newSP(1) ~= ALL_RESULT(i,26) || newSP(2) ~= ALL_RESULT(i,27) || newSP(3) ~= ALL_RESULT(i,29) || ... 
%             newSP(4) ~= ALL_RESULT(i,30) || newSP(5) ~= ALL_RESULT(i,32) || newSP(6) ~= ALL_RESULT(i,33)
%         disp('Re-running H-k-Vp stacking')
%         SP(1) = newSP(1); SP(2) = newSP(2); SP(4) = newSP(3); SP(5) = newSP(4); SP(7) = newSP(5); SP(8) = newSP(6);
%         [tern_info, HkVp_RESULT, HkVp_good, Hk_good, std_err_HkVp, std_err_Hk] = HkVp(stack_ACCZ,stack_ACCR,stack_RFs,SP,rayp,zonly,AssVp,wtog,ptog,w,COL,fname2,fname21,fname3);
%         ALL_RESULT(i,:) = HkVp_RESULT;
%     end
    
    
end

%% WRITE OUT FINAL RESULTS
% THIS DOESN'T WRITE OUT TEXTFILE PROPERLY FOR SOME REASON. HAVE TO WRITE
% OUT SEPARATELY
% fidout = fopen(strcat(str1,'.txt'),'w');
% for i =1:size(ALL_RESULT,1)
%     fprintf(fidout,'%s %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3 %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.0f\n', ...
%         STALST{i},ALL_RESULT(i,1),ALL_RESULT(i,2),ALL_RESULT(i,3),ALL_RESULT(i,4),ALL_RESULT(i,5),ALL_RESULT(i,6),ALL_RESULT(i,7),ALL_RESULT(i,8),ALL_RESULT(i,9),ALL_RESULT(i,10),ALL_RESULT(i,11),ALL_RESULT(i,12),ALL_RESULT(i,13),ALL_RESULT(i,14),ALL_RESULT(i,15),ALL_RESULT(i,16),ALL_RESULT(i,17),ALL_RESULT(i,18),ALL_RESULT(i,19),ALL_RESULT(i,20),ALL_RESULT(i,21),ALL_RESULT(i,22),ALL_RESULT(i,23),ALL_RESULT(i,24),ALL_RESULT(i,25),ALL_RESULT(i,26),ALL_RESULT(i,27),ALL_RESULT(i,28),ALL_RESULT(i,29),ALL_RESULT(i,30),ALL_RESULT(i,31),ALL_RESULT(i,32),ALL_RESULT(i,33),ALL_RESULT(i,34),ALL_RESULT(i,35));
% end
% fclose(fidout);
fidout = fopen(strcat(str1,'_STAS.txt'),'w');
for i =1:size(ALL_RESULT,1)
    fprintf(fidout,'%s %f %f\n', STALST{i},STALOCS(i,1),STALOCS(i,2));
end
fclose(fidout);
dlmwrite(strcat(str1,'.txt'),ALL_RESULT,'delimiter',' ');
