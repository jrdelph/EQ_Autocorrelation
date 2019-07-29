%
% Cycle through stations to verify results and compare with other models
%
% UPDATE: 7/22/19
% This script pairs with Solver_hkVp_hk_Interactive.m
%   Previous version (StackPicking_error_contain.m) was standalone and part
%   of a multistep process. This script cuts down on the number of scripts
%   necessary to reach a final result.

close all; clear all


mainFold = pwd;
MATFILE = sprintf('ALL_RESULTS_inter_stderrfix.mat') % MATLAB file for storing results
load(MATFILE);
datafolder = sprintf('../RESULTS_INTER_0.68'); % Folder of *_terninfo.mat and *_results.mat files

tmp = strsplit(MATFILE,'.');
%str1 = strcat(tmp(1),'.',tmp(2),'.',tmp(3),'.',tmp(4));
str1 = strcat(tmp(1));

% Parameters in case HkVp stacking is necessary or solution is changed
toggle = 0; %Toggle for analysis and plotting of Linear (0) or PW stacked AC (1)
plotme = 0; % 1 if you want to print out PDF of waveforms (TAKES LONG TIME)
w = [0.4 0.2 0.4] % weighting in order for Ps, 3P1S, and Pmp
wtog = 0; % normalize ACC weighting by number of events in bin...
          % 0 if unweighted, 1 if weighted by # traces in bin
ptog = 0; % If ptog = 1, amplitudes of opposite sign predicted along moveout 
          % curve are not summed for stack. If 0, all amplitudes along moveout
          % are stacked (can lead to destructive amplitude summing)
AssVp = 6.2;
zonly = 1; % 1: only use Z autocorrelation amplitude data, 2: use Z and R
rayp = [0.0402, 0.0446, 0.0490, 0.0534, 0.0578, 0.0622, 0.0666, 0.0710, 0.0754, 0.0798];


lon = STALOCS(:,1);
lat = STALOCS(:,2);

%% LOAD MODELS FOR COMPARISON

% Crust1.0
CRST1 = load('crst1.0.xyz');
[X1 Y1 Z1] = xyz2grid(CRST1(:,1),CRST1(:,2),CRST1(:,3));
CRST1VP = load('crst1.0Vp.xyz');
VP = xyz2grid(CRST1VP(:,1),CRST1VP(:,2),CRST1VP(:,3));

% Schmandt H and Vs (Vp scaled by Brocher relationship)
SCH_H = load('H_Schmandt.xyz');
[X2 Y2 Z2] = xyz2grid(SCH_H(:,1),SCH_H(:,2),SCH_H(:,3));
%SCH_H = SCH_H((~isnan(SCH_H(:,3))),:);
SCH_VS = load('VsSchmandt.xyz');
[X3 Y3 VS] = xyz2grid(SCH_VS(:,1),SCH_VS(:,2),SCH_VS(:,3));
    

%% Make map of stations with data 
%  STILL NEED TO ADD COMPARISON XY PLOTS
minlat = min(lat)-1;
maxlat = max(lat)+1;
minlon = min(lon)-1;
maxlon = max(lon)+1;

% INTERACTIVE FIGURE
figstuff = get(0,'ScreenSize');
h_all = figure('Visible','on','Position',[ figstuff(3)/2 figstuff(4)*0.5 figstuff(3) figstuff(4)]);
ylim([-0.7 1.2]);
load coastlines

% % Thickness
% ax1 = subplot(3,3,[1 2]);
% %imagesc(X1,Y1,Z1)
% h = pcolor(X1,Y1,Z1);
% set(h,'EdgeColor','none');
% hold on
% plot(coastlon,coastlat,'k');
% axis equal
% axis([minlon maxlon minlat maxlat]);
% h2 = scatter(lon(RES(:,12)>1),lat(RES(:,12)>1),50,RES(RES(:,12)>1,2),'filled','MarkerEdgeColor',[ 1 1 1 ],'LineWidth',2);
% h2 = scatter(lon(RES(:,12)==1),lat(RES(:,12)==1),50,RES(RES(:,12)==1,2),'filled','MarkerEdgeColor',[ 0 0 0 ]);
% %h2 = scatter(lon,lat,50,RES(:,2),'filled','MarkerEdgeColor',[ 0 0 0 ]);
% h2.CDataSource = 'RES(:,2)';
% colorbar
% colormap(ax1,jet)
% caxis([10 50])

% Define lats, lons of stations
X21dat = lon(ALL_RESULT(:,35)>1);
X22dat = lon(ALL_RESULT(:,35)==1);
Y21dat = lat(ALL_RESULT(:,35)>1);
Y22dat = lat(ALL_RESULT(:,35)==1);

% Vp map
ax2 = subplot(3,3,[4 5]);
%imagesc(X1,Y1,Z1)
h = pcolor(X1,Y1,VP);
set(h,'EdgeColor','none');
hold on
plot(coastlon,coastlat,'k');
axis equal
axis([minlon maxlon minlat maxlat]);
VP31dat = ALL_RESULT(ALL_RESULT(:,35)>1,6);
VP32dat = ALL_RESULT(ALL_RESULT(:,35)==1,6);
h31 = scatter(X21dat,Y21dat,50,VP31dat,'filled','MarkerEdgeColor',[ 1 1 1 ],'LineWidth',1);
h32 = scatter(X22dat,Y22dat,50,VP32dat,'filled','MarkerEdgeColor',[ 0 0 0 ]);
h31.XDataSource = 'X21dat';
h32.XDataSource = 'X22dat';
h31.YDataSource = 'Y21dat';
h32.YDataSource = 'Y22dat';
h31.CDataSource = 'VP31dat';
h32.CDataSource = 'VP32dat';
colorbar
colormap(ax2,flipud(jet))
caxis([5.6 7.2])

% Vs map
Vs = ALL_RESULT(:,6)./ALL_RESULT(:,11);
ax3 = subplot(3,3,[7 8]);
%imagesc(X1,Y1,Z1)
h = pcolor(X1,Y1,VS);
set(h,'EdgeColor','none');
hold on
plot(coastlon,coastlat,'k');
axis equal
axis([minlon maxlon minlat maxlat]);
VS41dat = Vs(ALL_RESULT(:,35)>1);
VS42dat = Vs(ALL_RESULT(:,35)==1);
h41 = scatter(X21dat,Y21dat,50,VS41dat,'filled','MarkerEdgeColor',[ 1 1 1 ],'LineWidth',1);
h42 = scatter(X22dat,Y22dat,50,VS42dat,'filled','MarkerEdgeColor',[ 0 0 0 ]);
h41.XDataSource = 'X21dat';
h42.XDataSource = 'X22dat';
h41.YDataSource = 'Y21dat';
h42.YDataSource = 'Y22dat';
h41.CDataSource = 'VS41dat';
h42.CDataSource = 'VS42dat';
colorbar
colormap(ax3,flipud(jet))
caxis([3 4.2])

% PLOT DIFF CHARTS
for i = 1:size(STALOCS,1)
    % COMPARE TO CRUST 1.0
    r_crst = sqrt((STALOCS(i,1) - CRST1(:,1)).^2 + (STALOCS(i,2) - CRST1(:,2)).^2);
    SEL2 = find(r_crst==min(r_crst)); 
    C1_H_all(i) = CRST1(SEL2,3);
    r_crstvp = sqrt((STALOCS(i,1) - CRST1VP(:,1)).^2 + (STALOCS(i,2) - CRST1VP(:,2)).^2);
    SEL3 = find(r_crstvp==min(r_crstvp)); 
    C1_Vp_all(i) = CRST1VP(SEL3,3);
    % COMPARE TO SCHMANDT
    r_SCH = sqrt((STALOCS(i,1) - SCH_VS(:,1)).^2 + (STALOCS(i,2) - SCH_VS(:,2)).^2);
    SEL4 = find(r_SCH==min(r_SCH)); % find closest station to click
    S_Vs_all(i) = SCH_VS(SEL4,3);
    S_Vp_all(i) = 0.9409 + 2.0947.*S_Vs_all(i) - 0.8206.*S_Vs_all(i).^2 + 0.2683.*S_Vs_all(i).^3 - 0.0251.*S_Vs_all(i).^4;
    rH_SCH = sqrt((STALOCS(i,1) - SCH_H(:,1)).^2 + (STALOCS(i,2) - SCH_H(:,2)).^2);
    SEL5 = find(rH_SCH==min(rH_SCH)); % find closest station to click
    S_H_all(i) = SCH_H(SEL5,3);
    %SCHLON(i) = SCH_H(SEL5,1);
    %SCHLAT(i) = SCH_H(SEL5,2);
end

%% NOT SURE WHAT TO DO HERE
%Vs_err = RES(:,6)./RES(:,11).*sqrt((RES(:,6)./RES(:,5)).^2 + (RES(:,9)./RES(:,8)).^2);

subplot(3,3,3)
%scatter(C1_H,RES(:,2),30,'k','filled')
h5 = errorbar(C1_H_all,ALL_RESULT(:,1),ALL_RESULT(:,1)-ALL_RESULT(:,3),ALL_RESULT(:,4)-ALL_RESULT(:,1),'k.','MarkerSize',10);
h5.YDataSource = 'ALL_RESULT(:,1)';
h5.UDataSource = 'ALL_RESULT(:,4)-ALL_RESULT(:,1)';
h5.LDataSource = 'ALL_RESULT(:,1)-ALL_RESULT(:,3)';
grid on
hold on
RMSE_H1 = sqrt(mean((ALL_RESULT(:,1)-C1_H_all').^2));
h51t = text(22,58,sprintf('RMSdiff = %0.2f',RMSE_H1));
RMSE_H2 = sqrt(mean((ALL_RESULT(:,1)-S_H_all').^2));
h52t = text(22,53,sprintf('RMSdiff = %0.2f',RMSE_H2));
xlabel('Crust1.0 H');
ylabel('H-k-Vp H');
axis([20 60 20 60])


% TEST = [ STALOCS(:,1) STALOCS(:,2) RES(:,2) ]
% TEST2 = [ SCHLON' SCHLAT' S_H_all' ]
% plot(SCH_H(:,1),SCH_H(:,2),'.')
% hold on
% scatter(TEST(:,1),TEST(:,2),'k','filled')
% scatter(TEST2(:,1),TEST2(:,2),'r','filled')


subplot(3,3,6)
%scatter(RES(:,5),C1_Vp,30,'k','filled')
h6 = errorbar(C1_Vp_all,ALL_RESULT(:,6),ALL_RESULT(:,6)-ALL_RESULT(:,8),ALL_RESULT(:,9)-ALL_RESULT(:,6),'k.','MarkerSize',10);
h6.YDataSource = 'ALL_RESULT(:,6)';
h6.LDataSource = 'ALL_RESULT(:,6)-ALL_RESULT(:,8)';
h6.UDataSource = 'ALL_RESULT(:,9)-ALL_RESULT(:,6)';
grid on
RMSE_Vp = sqrt(mean((ALL_RESULT(:,6)-C1_Vp_all').^2));
h6t = text(5.55,5.45,sprintf('RMSdiff = %0.2f',RMSE_Vp));
xlabel('Crust1.0 Vp');
ylabel('H-k-Vp Vp');
axis([5.5 7.2 5.25 7.4])

subplot(3,3,9)
%scatter(RES(:,5)./RES(:,8),S_Vs,30,'k','filled')
%h7 = errorbar(S_Vs_all,Vs,Vs_err,'k.','MarkerSize',10);
h7 = scatter(S_Vs_all,Vs,20,'ok','filled');
h7.YDataSource = 'ALL_RESULT(:,6)./ALL_RESULT(:,11)';
%h7.LDataSource = 'Vs_err';
%h7.UDataSource = 'Vs_err';
grid on
RMSE_Vs = sqrt(mean((ALL_RESULT(:,6)./ALL_RESULT(:,11)-S_Vs_all').^2));
h7t = text(3.02,4.8,sprintf('RMSdiff = %0.2f',RMSE_Vs));
xlabel('S2015 Vs');
ylabel('H-k-Vp Vs');
axis([3 4.0 2.5 5])

% Thickness
ax1 = subplot(3,3,[1 2]);
%imagesc(X1,Y1,Z1)
h = pcolor(X1,Y1,Z1);
set(h,'EdgeColor','none');
hold on
plot(coastlon,coastlat,'k');
axis equal
axis([minlon maxlon minlat maxlat]);
H21dat = ALL_RESULT(ALL_RESULT(:,35)>1,1);
H22dat = ALL_RESULT(ALL_RESULT(:,35)==1,1);
h21 = scatter(X21dat,Y21dat,50,H21dat,'filled','MarkerEdgeColor',[ 1 1 1 ],'LineWidth',1);
h22 = scatter(X22dat,Y22dat,50,H22dat,'filled','MarkerEdgeColor',[ 0 0 0 ]);
h21.XDataSource = 'X21dat';
h22.XDataSource = 'X22dat';
h21.YDataSource = 'Y21dat';
h22.YDataSource = 'Y22dat';
h21.CDataSource = 'H21dat';
h22.CDataSource = 'H22dat';
colorbar
colormap(ax1,jet)
caxis([10 50])


%%
loop1 = true;

while loop1

    % Select point
    [x1,y1,button1] = ginput(1)
    if button1 == 1
        r = sqrt((x1 - lon).^2 + (y1 - lat).^2);
        SEL = find(r==min(r)); % find closest station to click
        STA_SEL = STALST{SEL}
        STALON = lon(SEL);
        STALAT = lat(SEL);
        
        r_crst = sqrt((STALON - CRST1(:,1)).^2 + (STALAT - CRST1(:,2)).^2);
        SEL2 = find(r_crst==min(r_crst)); % find closest station to click
        C1_H = CRST1(SEL2,3);
        C1_Vp = CRST1VP(SEL2,3);
   
        r_SCH = sqrt((STALON - SCH_VS(:,1)).^2 + (STALAT - SCH_VS(:,2)).^2);
        SEL4 = find(r_SCH==min(r_SCH)); % find closest station to click
        S_Vs = SCH_VS(SEL4,3);
        S_Vp = 0.9409 + 2.0947.*S_Vs - 0.8206.*S_Vs.^2 + 0.2683.*S_Vs.^3 - 0.0251.*S_Vs.^4;
        r_SCH = sqrt((STALON - SCH_H(:,1)).^2 + (STALAT - SCH_H(:,2)).^2);
        SEL4 = find(r_SCH==min(r_SCH)); % find closest station to click
        S_H = SCH_H(SEL4,3);

        % LOAD DATA FROM CLICK
        f = strcat('./',datafolder,'/',STA_SEL,'_results.mat');
        fname2 = strcat('./',datafolder,'/',STA_SEL,'_terninfo.txt');
        fname21 = strcat('./',datafolder,'/',STA_SEL,'_terninfo.mat');
        fname3 = strcat('./',datafolder,'/',STA_SEL,'_BHinfo.txt');
        fname4 = strcat('./',datafolder,'/',STA_SEL,'_clustsolns.txt');

        load(f);
        load(fname21);
        idx = strmatch(STA_SEL,STALST);

        % PLOT  RESULTS AND CHECK SOLUTIONS
        % This function allows you to change parameters for HkVp/Hk
        % stacking if you want, but doesn't run HkVp stacking itself.
        % Cluster analysis, however, is run in this script.
        GOOD = 0;
        if toggle == 0
            COL = 5;
            fname = strcat('./',datafolder,'/',char(STA_SEL),'_HVpk_LIN.pdf');
            scale = 0.005;
        elseif toggle == 1
            COL = 6;
            fname = strcat('./',datafolder,'/',char(STA_SEL),'_HVpk_PW.pdf');
            scale = 0.008;
        else
            disp('Stacking not specified. Using linear stacking')
            COL = 5;
            fname = strcat('./',datafolder,'/',char(STA_SEL),'_HVpk_LIN.pdf');
            scale = 0.005;
        end
        % Until good solution is declared
        while GOOD == 0
            tic
            % Can delete first HkVp function in future version as it is now saved in
            % *_terninfo.mat from previous script
            w = ALL_RESULT(idx,36:38);
            [tern_info, HkVp_RESULT, HkVp_good, Hk_good, std_err_HkVp, std_err_Hk] = HkVp(stack_ACCZ,stack_ACCR,stack_RFs,ALL_RESULT(idx,26:34),rayp,zonly,AssVp,wtog,ptog,w,COL,fname2,fname21,fname3);
            ALL_RESULT(idx,1)=HkVp_RESULT(1);ALL_RESULT(idx,2)=HkVp_RESULT(2);ALL_RESULT(idx,3)=HkVp_RESULT(3);ALL_RESULT(idx,4)=HkVp_RESULT(4);ALL_RESULT(idx,5)=HkVp_RESULT(5);ALL_RESULT(idx,6)=HkVp_RESULT(6);
            ALL_RESULT(idx,7)=HkVp_RESULT(7);ALL_RESULT(idx,8)=HkVp_RESULT(8);ALL_RESULT(idx,9)=HkVp_RESULT(9);ALL_RESULT(idx,10)=HkVp_RESULT(10);ALL_RESULT(idx,11)=HkVp_RESULT(11);ALL_RESULT(idx,12)=HkVp_RESULT(12);
            ALL_RESULT(idx,13)=HkVp_RESULT(13);ALL_RESULT(idx,14)=HkVp_RESULT(14);ALL_RESULT(idx,15)=HkVp_RESULT(15);ALL_RESULT(idx,16)=HkVp_RESULT(16);ALL_RESULT(idx,17)=HkVp_RESULT(17);ALL_RESULT(idx,18)=HkVp_RESULT(18);
            ALL_RESULT(idx,19)=HkVp_RESULT(19);ALL_RESULT(idx,20)=HkVp_RESULT(20);ALL_RESULT(idx,21)=HkVp_RESULT(21);ALL_RESULT(idx,22)=HkVp_RESULT(22);ALL_RESULT(idx,23)=HkVp_RESULT(23);ALL_RESULT(idx,24)=HkVp_RESULT(24);    
            ALL_RESULT(idx,25)=HkVp_RESULT(25);ALL_RESULT(idx,26)=HkVp_RESULT(26);ALL_RESULT(idx,27)=HkVp_RESULT(27);ALL_RESULT(idx,28)=HkVp_RESULT(28);ALL_RESULT(idx,29)=HkVp_RESULT(29);ALL_RESULT(idx,30)=HkVp_RESULT(30);
            ALL_RESULT(idx,31)=HkVp_RESULT(31);ALL_RESULT(idx,32)=HkVp_RESULT(32);ALL_RESULT(idx,33)=HkVp_RESULT(33);ALL_RESULT(idx,34)=HkVp_RESULT(34);ALL_RESULT(idx,36)=HkVp_RESULT(36);ALL_RESULT(idx,37)=HkVp_RESULT(37);
            ALL_RESULT(idx,38)=HkVp_RESULT(38);
            
            [GOOD, newSOLN] =plot_HkVp(fname,fname2,fname4,stack_ACCZ,stack_ACCR,stack_RFs,ALL_RESULT(idx,:),std_err_HkVp,std_err_Hk,HkVp_good,Hk_good,COL,scale,plotme);
            %  This "if" resets n_clust to 1 if boundaries change for hkVp
            %  stacking. Retains grouping if not.
            if sum(HkVp_RESULT(26:34) ~= newSOLN(26:34)) > 0 || sum(HkVp_RESULT(36:38)~= newSOLN(36:38))
                newSOLN(35) = 1;
            end
            ALL_RESULT(idx,:) = newSOLN;
        end
        save(char(strcat(str1,'_mappick.mat')),'STALST','STALOCS','ALL_RESULT');
    
        %[newsoln, GOOD, groups] = plot_results_noOUT(stack_ACCZ,stack_ACCR,stack_RFs,1,pdf_ans,tern_info,[C1_H C1_Vp],[S_H S_Vp S_Vs],0.95);
        figure(1)

        % WRITE OUT SOLUTION TO FILE
        %staidx = find(strcmp(STA_SEL,cellstr(STALST))==1);
        %ALL_RESULT(idx,:) = newSOLN; % Replace old solution
        % Recalc Vs_err
        Vs = ALL_RESULT(:,6)./ALL_RESULT(:,11);
        %Vs_err = ALL_RESULT(:,6)./ALL_RESULT(:,11).*sqrt((ALL_RESULT(:,6)./ALL_RESULT(:,5)).^2 + (ALL_RESULT(:,9)./ALL_RESULT(:,8)).^2);


        delete(h51t)
        delete(h52t)
        subplot(3,3,3)
        RMSE_H1 = sqrt(mean((ALL_RESULT(:,1)-C1_H_all').^2));
        RMSE_H2 = sqrt(mean((ALL_RESULT(:,1)-S_H_all').^2));
        h51t = text(22,55,sprintf('RMSdiff = %0.2f',RMSE_H1));
        h52t = text(22,53,sprintf('RMSdiff = %0.2f',RMSE_H2));
        % Get percentage of measurements that fall within error
        DIFF_C1 = abs((ALL_RESULT(:,1)-C1_H_all'));
        PERC_C1 = sum(DIFF_C1<5)/length(DIFF_C1);
        DIFF_S1 = abs((ALL_RESULT(:,1)-S_H_all'));
        PERC_S1 = sum(DIFF_S1<5)/length(DIFF_S1);
        
        RMSE_Vp = sqrt(mean((ALL_RESULT(:,6)-C1_Vp_all').^2));
        delete(h6t)
        subplot(3,3,6)
        h6t = text(5.55,5.5,sprintf('RMSdiff = %0.2f',RMSE_Vp));
        RMSE_Vs = sqrt(mean((Vs-S_Vs_all').^2));
        delete(h7t)
        subplot(3,3,9)
        h7t = text(3.02,4.8,sprintf('RMSdiff = %0.2f',RMSE_Vs));
        
        clear X21dat X22dat Y21dat Y22dat H21dat H22dat VP31dat VP32dat VS41dat VS42dat
        X21dat = lon(ALL_RESULT(:,35)>1);
        X22dat = lon(ALL_RESULT(:,35)==1);
        Y21dat = lat(ALL_RESULT(:,35)>1);
        Y22dat = lat(ALL_RESULT(:,35)==1);
        H21dat = ALL_RESULT(ALL_RESULT(:,35)>1,1);
        H22dat = ALL_RESULT(ALL_RESULT(:,35)==1,1);  
        VP31dat = ALL_RESULT(ALL_RESULT(:,35)>1,6);
        VP32dat = ALL_RESULT(ALL_RESULT(:,35)==1,6);
        VS41dat = Vs(ALL_RESULT(:,35)>1);
        VS42dat = Vs(ALL_RESULT(:,35)==1);
        
        % Update result comparison
        figure(1)
        refreshdata
        h_all.CurrentAxes = ax1;
        
        %save('Solutions.mat','STALST','STALOCS','ALL_RESULT');

        clear stack_ACCZ stack_ACCR stack_RFs pdf_ans tern_info
        elseif button1 ==3
            outfile = cell2mat(inputdlg({'Enter resultfile name'},'Export',1,{'OUT_results.txt'}));
            fidout = fopen(outfile,'w');
            fprintf(fidout,'STNM LON LAT H_best H_min H_25 H_75 H_max Vp_best Vp_min Vp_25 Vp_75 Vp_max VpVs_best VpVs_min VpVs_25 VpVs_75 VpVs_max H_Hkbest H_Hkmin H_Hk25 H_Hk75 H_Hkmax VpVs_Hkbest VpVs_Hkmin VpVs_Hk25 VpVs_Hk75 VpVs_Hkmax Hsearch_min Hsearch_max Hsearch_step Vpsearch_min Vpsearch_max Vpsearch_step VpVssearch_min VpVssearch_max VpVssearch_step n_clust weight_Ps weight_2P1S weight_Pmp\n');          
            for i = 1:size(ALL_RESULT,1)
                fprintf(fidout,'%s %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.0f %0.2f %0.2f %0.2f\n',STALST{i}, STALOCS(i,1),STALOCS(i,2),...
                    ALL_RESULT(i,1),ALL_RESULT(i,2),ALL_RESULT(i,3),ALL_RESULT(i,4),ALL_RESULT(i,5),ALL_RESULT(i,6),ALL_RESULT(i,7),ALL_RESULT(i,8),ALL_RESULT(i,9),ALL_RESULT(i,10),ALL_RESULT(i,11),ALL_RESULT(i,12),ALL_RESULT(i,13),ALL_RESULT(i,14),...
                    ALL_RESULT(i,15),ALL_RESULT(i,16),ALL_RESULT(i,17),ALL_RESULT(i,18),ALL_RESULT(i,19),ALL_RESULT(i,20),ALL_RESULT(i,21),ALL_RESULT(i,22),ALL_RESULT(i,23),ALL_RESULT(i,24),ALL_RESULT(i,25),ALL_RESULT(i,26),ALL_RESULT(i,27),ALL_RESULT(i,28), ...
                    ALL_RESULT(i,29),ALL_RESULT(i,30),ALL_RESULT(i,31),ALL_RESULT(i,32),ALL_RESULT(i,33),ALL_RESULT(i,34),ALL_RESULT(i,35),ALL_RESULT(i,36),ALL_RESULT(i,37),ALL_RESULT(i,38));
            end
            % OUTPUT VS VSERR
            dlmwrite('H_comp.txt',[ALL_RESULT(:,1) ALL_RESULT(:,3) ALL_RESULT(:,4) C1_H_all'],'delimiter',' ');
            dlmwrite('Vp_comp.txt',[ALL_RESULT(:,6) ALL_RESULT(:,8) ALL_RESULT(:,9) C1_Vp_all'],'delimiter',' ');
            dlmwrite('Vs_comp.txt',[Vs S_Vs_all'],'delimiter',' ');
            
            loop1=false;
        end 


        
    end





