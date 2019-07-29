function [GOOD,HkVp_RESULT] = plot_HkVp(fname,fname2,fname4,stack_ACCZ,stack_ACCR,stack_RFs,HkVp_RESULT,std_err_HkVp,std_err_Hk,HkVp_good,Hk_good,COL,scale,plotme)
% This function creates a plot of all waveforms as well as travel times from 
% H-k and H-k-Vp solutions.
%
% fname is output pdf name (if plotme = 1)
% stack_* are binned and stacked waveforms from EQ_AUTOCORR
% HkVp_good and Hk_good are solutions within standard error of maximum
% COL define whether linear or pw stacks are plotted (5 = linear, 6 = pwstack)

n_clust = HkVp_RESULT(35);

AssVp = 6.2;

rayp_pred = 0.00:0.002:0.10;

H_Hkbest = HkVp_RESULT(16);
H_Hkq = [ HkVp_RESULT(18) HkVp_RESULT(19) ];

VpVs_Hkbest = HkVp_RESULT(21);
VpVs_Hkq = [ HkVp_RESULT(23) HkVp_RESULT(24) ];

minH = HkVp_RESULT(26);    
maxH = HkVp_RESULT(27);  
stepH = HkVp_RESULT(28);  
minVp = HkVp_RESULT(29);    
maxVp= HkVp_RESULT(30); 
stepVp = HkVp_RESULT(31);  
minVpVs = HkVp_RESULT(32);    
maxVpVs = HkVp_RESULT(33);
stepVpVs = HkVp_RESULT(34);  
w = [HkVp_RESULT(36) HkVp_RESULT(37) HkVp_RESULT(38)];

H_good = HkVp_good(:,1);
Vp_good = HkVp_good(:,2);
VpVs_good = HkVp_good(:,3);
Amps_good = HkVp_good(:,4);

X = [(H_good - minH)./(maxH-minH), ...
    (Vp_good - minVp)./(maxVp-minVp), ...
    (VpVs_good - minVpVs)./(maxVpVs-minVpVs), ...
    (Amps_good - std_err_HkVp)./(max(Amps_good)-std_err_HkVp)];

% 
Tpmp_good = zeros(length(H_good),length(rayp_pred));
Tps_good = zeros(length(H_good),length(rayp_pred));
Tppps_good = zeros(length(H_good),length(rayp_pred));
Tpsps_good = zeros(length(H_good),length(rayp_pred));

% Timing of all solutions with amplitude withing standard error of max
for ii = 1:length(H_good)
    Tpmp_good(ii,:) = 2*sqrt(1/Vp_good(ii).^2 - rayp_pred.^2)*H_good(ii);
    Tps_good(ii,:) = (sqrt(1/(Vp_good(ii)/VpVs_good(ii)).^2 - rayp_pred.^2) - sqrt(1/Vp_good(ii).^2 - rayp_pred.^2))*H_good(ii);
    Tppps_good(ii,:) = (sqrt(1/(Vp_good(ii)/VpVs_good(ii)).^2-rayp_pred.^2) + sqrt(1/Vp_good(ii).^2 - rayp_pred.^2))*H_good(ii);
    Tpsps_good(ii,:) = 2*(sqrt(1/(Vp_good(ii)/VpVs_good(ii)).^2-rayp_pred.^2))*H_good(ii);
end

% GET TIMES FOR H-k soln
Tpmp_Hk = 2*sqrt(1/AssVp.^2 - rayp_pred.^2)*H_Hkbest;
Tps_Hk = (sqrt(1/(AssVp/VpVs_Hkbest).^2 - rayp_pred.^2) - sqrt(1/AssVp.^2 - rayp_pred.^2))*H_Hkbest;
Tppps_Hk = (sqrt(1/(AssVp/VpVs_Hkbest).^2-rayp_pred.^2) + sqrt(1/AssVp.^2 - rayp_pred.^2))*H_Hkbest;
Tpsps_Hk = 2*(sqrt(1/(AssVp/VpVs_Hkbest).^2-rayp_pred.^2))*H_Hkbest;
    
%% Parallel plots and histograms for single solutions
if n_clust == 1
    % Get best solution
    H_best = HkVp_RESULT(1);
    H_q = [ HkVp_RESULT(3) HkVp_RESULT(4) ];

    Vp_best = HkVp_RESULT(6);
    Vp_q = [ HkVp_RESULT(8) HkVp_RESULT(9) ];

    VpVs_best = HkVp_RESULT(11);
    VpVs_q = [ HkVp_RESULT(13) HkVp_RESULT(14) ];

    % GET TIMES FOR H-k-Vp soln
    Tpmp_pred = 2*sqrt(1/Vp_best.^2 - rayp_pred.^2)*H_best;
    Tps_pred = (sqrt(1/(Vp_best/VpVs_best).^2 - rayp_pred.^2) - sqrt(1/Vp_best.^2 - rayp_pred.^2))*H_best;
    Tppps_pred = (sqrt(1/(Vp_best/VpVs_best).^2-rayp_pred.^2) + sqrt(1/Vp_best.^2 - rayp_pred.^2))*H_best;
    Tpsps_pred = 2*(sqrt(1/(Vp_best/VpVs_best).^2-rayp_pred.^2))*H_best;

    % Make figures

    % PARALLEL PLOT
    h = figure(4);
    %h.Renderer='Painters';
    set(h,'pos',[318 598 934 528]);

    minH = HkVp_RESULT(26);    
    maxH = HkVp_RESULT(27);  
    stepH = HkVp_RESULT(28);  
    minVp = HkVp_RESULT(29);    
    maxVp= HkVp_RESULT(30); 
    stepVp = HkVp_RESULT(31);  
    minVpVs = HkVp_RESULT(32);    
    maxVpVs = HkVp_RESULT(33);
    stepVpVs = HkVp_RESULT(34);  

    X = [(H_good - minH)./(maxH-minH), ...
        (Vp_good - minVp)./(maxVp-minVp), ...
        (VpVs_good - minVpVs)./(maxVpVs-minVpVs), ...
        (Amps_good - std_err_HkVp)./(max(Amps_good)-std_err_HkVp)];

    idx = ones(size(X,1),1);


    subplot(1,6,[1 2 3])
    h1 = parallelcoords(X,'Group',idx,'Labels',{'H','Vp','VpVs','Amp'});
    grid on
    ylim([0 1])
    
    subplot(1,6,4)
    hold on
    h_temp = histogram(H_good,minH:stepH:maxH,'orientation','horizontal');
    plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[HkVp_RESULT(2),HkVp_RESULT(5)],'k','LineWidth',2)
    if HkVp_RESULT(4)-HkVp_RESULT(3)<1e-8
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(3) (max(h_temp.Values)*1.1)/10 0.5],'FaceColor','black')
    else
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(3) (max(h_temp.Values)*1.1)/10 HkVp_RESULT(4)-HkVp_RESULT(3)],'FaceColor','black')
    end
    plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[HkVp_RESULT(1),HkVp_RESULT(1)],'y','LineWidth',3)
    axis([ 0 max(h_temp.Values)*1.5 minH maxH])
    
    subplot(1,6,5)
    hold on
    h_temp = histogram(Vp_good,minVp:stepVp:maxVp,'orientation','horizontal');
    plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[HkVp_RESULT(7),HkVp_RESULT(10)],'k','LineWidth',2)
    % In rare case that range of solutions is 0, this prevents plotting error
    if HkVp_RESULT(9)-HkVp_RESULT(8)<1e-8
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(8) (max(h_temp.Values)*1.1)/10 0.5],'FaceColor','black')
    else
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(8) (max(h_temp.Values)*1.1)/10 HkVp_RESULT(9)-HkVp_RESULT(8)],'FaceColor','black')
    end
    plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[HkVp_RESULT(6),HkVp_RESULT(6)],'y','LineWidth',3)
    axis([ 0 max(h_temp.Values)*1.5 minVp maxVp])
    
    subplot(1,6,6)
    hold on
    h_temp = histogram(VpVs_good,minVpVs:stepVpVs:maxVpVs,'orientation','horizontal');
    plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[HkVp_RESULT(12),HkVp_RESULT(15)],'k','LineWidth',2)
    % In rare case that range of solutions is 0, this prevents plotting error
    if HkVp_RESULT(14)-HkVp_RESULT(13)<1e-8
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(13) (max(h_temp.Values)*1.1)/10 0.5],'FaceColor','black')
    else
        rectangle('Position',[max(h_temp.Values)*1.1 HkVp_RESULT(13) (max(h_temp.Values)*1.1)/10 HkVp_RESULT(14)-HkVp_RESULT(13)],'FaceColor','black')
    end
    plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[HkVp_RESULT(11),HkVp_RESULT(11)],'y','LineWidth',3)   
    axis([ 0 max(h_temp.Values)*1.5 minVpVs maxVpVs])
    
end
        
       
%% Parallel plots and histograms for multiple solutions
%  Also interactive selection for best pick

if n_clust > 1
    
    % RUN CLUSTER ANALYSIS
    % Gaussian Mixture modeling
    cc = lines(n_clust);
    clear idx
    gm = fitgmdist(X(:,1:3),n_clust,'RegularizationValue',0.001)
    idx = cluster(gm,X(:,1:3));

    % PARALLEL PLOT
    h = figure(4);
    %h.Renderer='Painters';
    set(h,'pos',[318 598 934 528]);
    
    subplot(1,6,[1 2 3])
    h1 = parallelcoords(X,'Group',idx,'Labels',{'H','Vp','VpVs','Amp'});
    grid on
    ylim([0 1])
    
    subplot(1,6,4)
    hold on
    for ii = 1:n_clust
        h_temp = histogram(H_good(idx==ii),minH:stepH:maxH,'orientation','horizontal','FaceColor',cc(ii,:));
        H_q(ii,:) = quantile(H_good(idx==ii),[.159 .841]); % Find middle 68% of distribution for each solution
        plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[min(H_good(idx==ii)),max(H_good(idx==ii))],'color',cc(ii,:),'LineWidth',2)
        rectangle('Position',[max(h_temp.Values)*1.1 H_q(ii,1) (max(h_temp.Values)*1.1)/10 H_q(ii,2)-H_q(ii,1)],'FaceColor',cc(ii,:))
        plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[H_good(Amps_good==max(Amps_good(idx==ii))),H_good(Amps_good==max(Amps_good(idx==ii)))],'y','LineWidth',3) 
        HSOLNS(ii,:) = [H_good(Amps_good==max(Amps_good(idx==ii))) H_q(ii,1) H_q(ii,2) min(H_good(idx==ii)) max(H_good(idx==ii)) ];
    end
    axis([ 0 max(h_temp.Values)*1.5 minH maxH])

    subplot(1,6,5)
    hold on
    for ii = 1:n_clust
        h_temp = histogram(Vp_good(idx==ii),minVp:stepVp:maxVp,'orientation','horizontal','FaceColor',cc(ii,:));
        Vp_q(ii,:) = quantile(Vp_good(idx==ii),[.159 .841]); % Find middle 68% of distribution for each solution
        plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[min(Vp_good(idx==ii)),max(Vp_good(idx==ii))],'color',cc(ii,:),'LineWidth',2)
        rectangle('Position',[max(h_temp.Values)*1.1 Vp_q(ii,1) (max(h_temp.Values)*1.1)/10 Vp_q(ii,2)-Vp_q(ii,1)],'FaceColor',cc(ii,:))
        plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[Vp_good(Amps_good==max(Amps_good(idx==ii))),Vp_good(Amps_good==max(Amps_good(idx==ii)))],'y','LineWidth',3) 
        VPSOLNS(ii,:) = [Vp_good(Amps_good==max(Amps_good(idx==ii))) Vp_q(ii,1) Vp_q(ii,2) min(Vp_good(idx==ii)) max(Vp_good(idx==ii))];    
    end
    axis([ 0 max(h_temp.Values)*1.5 minVp maxVp])
    
    subplot(1,6,6)
    hold on
    for ii = 1:n_clust
        h_temp = histogram(VpVs_good(idx==ii),minVpVs:stepVpVs:maxVpVs,'orientation','horizontal','FaceColor',cc(ii,:));
        VpVs_q(ii,:) = quantile(VpVs_good(idx==ii),[.159 .841]); % Find middle 68% of distribution for each solution
        plot([max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)/2],[min(VpVs_good(idx==ii)),max(VpVs_good(idx==ii))],'color',cc(ii,:),'LineWidth',2)
        rectangle('Position',[max(h_temp.Values)*1.1 VpVs_q(ii,1) (max(h_temp.Values)*1.1)/10 VpVs_q(ii,2)-VpVs_q(ii,1)],'FaceColor',cc(ii,:))
        plot([max(h_temp.Values)*1.1,max(h_temp.Values)*1.1 + ((max(h_temp.Values)*1.1)/10)],[VpVs_good(Amps_good==max(Amps_good(idx==ii))),VpVs_good(Amps_good==max(Amps_good(idx==ii)))],'y','LineWidth',3) 
        VPVSSOLNS(ii,:) = [VpVs_good(Amps_good==max(Amps_good(idx==ii))) VpVs_q(ii,1) VpVs_q(ii,2) min(VpVs_good(idx==ii)) max(VpVs_good(idx==ii))];        
    end
    axis([ 0 max(h_temp.Values)*1.5 minVpVs maxVpVs])
    
    % Write out cluster solutions
    dlmwrite(fname4,[ [1:n_clust]' HSOLNS VPSOLNS VPVSSOLNS ],'delimiter',' ');
    
    h = figure(5);
    %h.Renderer='Painters';
    set(h,'pos',[200 300 934 528]);
    tshift = 0.02;
    tcent = 0.6;
    for ii = 1:n_clust
        if ii > 1
            tcent = tcent - 0.10;
        end
        text(0.05,tcent,sprintf('Solution %0.0f:',ii),'FontSize',16);
        text(0.2,tcent,sprintf('H = %0.1f',HSOLNS(ii,1)),'FontSize',16);
        text(0.28,tcent-tshift,sprintf('+ %0.1f',HSOLNS(ii,3)-HSOLNS(ii,1)),'FontSize',12);
        text(0.28,tcent+tshift,sprintf('- %0.1f',HSOLNS(ii,1)-HSOLNS(ii,2)),'FontSize',12);
        text(0.41,tcent,sprintf('Vp = %0.2f',VPSOLNS(ii,1)),'FontSize',16);
        text(0.49,tcent-tshift,sprintf('+ %0.2f',VPSOLNS(ii,3)-VPSOLNS(ii,1)),'FontSize',12);
        text(0.49,tcent+tshift,sprintf('- %0.2f',VPSOLNS(ii,1)-VPSOLNS(ii,2)),'FontSize',12);
        text(0.62,tcent,sprintf('VpVs = %0.3f',VPVSSOLNS(ii,1)),'FontSize',16);
        text(0.7,tcent-tshift,sprintf('+ %0.3f',VPVSSOLNS(ii,3)-VPVSSOLNS(ii,1)),'FontSize',12);
        text(0.7,tcent+tshift,sprintf('- %0.3f',VPVSSOLNS(ii,1)-VPVSSOLNS(ii,2)),'FontSize',12);        
    end

end


    
%% WAVEFORMS AND SOLUTIONS
h = figure(3);
h.Renderer='Painters';

%h = figure(3,'visible','off','pos',[500 500 1640 1000]);
%set(h,'pos',[500 500 656 400],'visible','off','PaperSize',[11 11]);
set(h,'pos',[1130 392 1042 757],'PaperSize',[15 15]);

t_width = 0.3;
t_height = 0.7;
h_width = t_width;
h_height = 0.15;

% VERTICAL AUTOCORRELATIONS
% HISTOGRAM 1
subplot(3,3,1,'Position',[0.02,0.75,h_width,h_height])
bar(cell2mat(stack_ACCZ(:,1)),cell2mat(stack_ACCZ(:,3)))
max_y = max(cell2mat(stack_ACCZ(:,3)))+.2*max(cell2mat(stack_ACCZ(:,3)));
title('Z AUTOCORR')
set(gca,'xticklabel',{});
axis([0.03 0.09 0 max_y])
for ii = 1:size(stack_ACCZ,1)
    text(cell2mat(stack_ACCZ(ii,1)),cell2mat(stack_ACCZ(ii,3))+.1*max_y, ...
        sprintf('%0.0f',cell2mat(stack_ACCZ(ii,3))),'VerticalAlignment','top','HorizontalAlignment','center')
end

subplot(3,3,[4 7],'Position',[0.02,0.05,t_width,t_height])
hold on
if n_clust > 1
    for ii = 1:n_clust
        plot(rayp_pred,Tpmp_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tppps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tpsps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
    end
    for ii = 1:size(stack_ACCZ,1)
        Time = [min(stack_ACCZ{ii,4})', stack_ACCZ{ii,4}', max(stack_ACCZ{ii,4})'];
        posAmps = [stack_ACCZ{ii,1}', max((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,1})', stack_ACCZ{ii,1}'];
        negAmps = [stack_ACCZ{ii,1}', min((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,1})', stack_ACCZ{ii,1}'];
        fill([posAmps,zeros(size(posAmps))+posAmps(1)],[Time,flip(Time)],'r','LineStyle','none')
        fill([negAmps,zeros(size(negAmps))+negAmps(1)],[Time,flip(Time)],'b','LineStyle','none')
        plot((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,4},'k')
    end
    for ii = 1:n_clust
        % GET TIMES FOR H-k-Vp solns
        Tpmp_pred = 2*sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2)*HSOLNS(ii,1);
        Tps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2 - rayp_pred.^2) - sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tppps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2) + sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tpsps_pred = 2*(sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2))*HSOLNS(ii,1);
        plot(rayp_pred,Tpmp_pred,'b','LineWidth',2)
        plot(rayp_pred,Tps_pred,'b','LineWidth',2)
        plot(rayp_pred,Tppps_pred,'b','LineWidth',2)
        plot(rayp_pred,Tpsps_pred,'r','LineWidth',2)
    end
else
    plot(rayp_pred,Tpmp_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tppps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tpsps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    for ii = 1:size(stack_ACCZ,1)
        Time = [min(stack_ACCZ{ii,4})', stack_ACCZ{ii,4}', max(stack_ACCZ{ii,4})'];
        posAmps = [stack_ACCZ{ii,1}', max((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,1})', stack_ACCZ{ii,1}'];
        negAmps = [stack_ACCZ{ii,1}', min((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,1})', stack_ACCZ{ii,1}'];
        fill([posAmps,zeros(size(posAmps))+posAmps(1)],[Time,flip(Time)],'r','LineStyle','none')
        fill([negAmps,zeros(size(negAmps))+negAmps(1)],[Time,flip(Time)],'b','LineStyle','none')
        plot((stack_ACCZ{ii,COL}./max(abs(stack_ACCZ{ii,COL}))).*scale+stack_ACCZ{ii,1},stack_ACCZ{ii,4},'k')
    end
    plot(rayp_pred,Tpmp_pred,'b','LineWidth',2)
    plot(rayp_pred,Tps_pred,'b','LineWidth',2)
    plot(rayp_pred,Tppps_pred,'b','LineWidth',2)
    plot(rayp_pred,Tpsps_pred,'r','LineWidth',2)
end
plot(rayp_pred,Tpmp_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tppps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tpsps_Hk,'k--','LineWidth',2)
axis([0.03 0.09 0 25])
axis ij
ylabel('Time (s)')
xlabel('Ray parameter (s/km)')

if n_clust == 1
    % TEXT OF ANSWERS
    tshift = 0.3;
    tcent = -7.2;
    text(0.031,tcent,sprintf('H-k-Vp Stacking:'));
    text(0.07,tcent,sprintf('H = %0.1f',H_best));
    text(0.082,tcent-tshift,sprintf('+ %0.1f',H_q(2)-H_best),'FontSize',8);
    text(0.082,tcent+tshift,sprintf('- %0.1f',H_best-H_q(1)),'FontSize',8);
    text(0.095,tcent,sprintf('Vp/Vs = %0.3f',VpVs_best));
    text(0.115,tcent-tshift,sprintf('+ %0.3f',VpVs_q(2)-VpVs_best),'FontSize',8);
    text(0.115,tcent+tshift,sprintf('- %0.3f',VpVs_best - VpVs_q(1)),'FontSize',8);                
    text(0.13,tcent,sprintf('Vp = %0.2f',Vp_best));
    text(0.143,tcent-tshift,sprintf('+ %0.2f',Vp_q(2)-Vp_best),'FontSize',8);
    text(0.143,tcent+tshift,sprintf('- %0.2f',Vp_best-Vp_q(1)),'FontSize',8);    
end

tshift = 0.3;
tcent = -8.4;
text(0.031,tcent,sprintf('H-k Stacking:'));
text(0.07,tcent,sprintf('H = %0.1f',H_Hkbest));
text(0.082,tcent-tshift,sprintf('+ %0.1f',H_Hkq(2)-H_Hkbest),'FontSize',8);
text(0.082,tcent+tshift,sprintf('- %0.1f',H_Hkbest-H_Hkq(1)),'FontSize',8);
text(0.095,tcent,sprintf('Vp/Vs = %0.3f',VpVs_Hkbest));
text(0.115,tcent-tshift,sprintf('+ %0.3f',VpVs_Hkq(2)-VpVs_Hkbest),'FontSize',8);
text(0.115,tcent+tshift,sprintf('- %0.3f',VpVs_Hkbest-VpVs_Hkq(1)),'FontSize',8);         
text(0.13,tcent,sprintf('Vp = %0.2f',AssVp));


% RADIAL 
subplot(3,3,2,'Position',[0.355,0.75,h_width,h_height])
bar(cell2mat(stack_ACCZ(:,1)),cell2mat(stack_ACCZ(:,3)))
max_y = max(cell2mat(stack_ACCZ(:,3)))+.2*max(cell2mat(stack_ACCZ(:,3)));
title('R AUTOCORR')
set(gca,'xticklabel',{});
axis([0.03 0.09 0 max_y])
for ii = 1:size(stack_ACCZ,1)
    text(cell2mat(stack_ACCZ(ii,1)),cell2mat(stack_ACCZ(ii,3))+.1*max_y, ...
        sprintf('%0.0f',cell2mat(stack_ACCZ(ii,3))),'VerticalAlignment','top','HorizontalAlignment','center')
end

subplot(3,3,[5 8],'Position',[0.355,0.05,t_width,t_height])
hold on
if n_clust > 1
    for ii = 1:n_clust
        plot(rayp_pred,Tpmp_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tppps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tpsps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
    end
    for ii = 1:size(stack_ACCR,1)
        Time = [min(stack_ACCR{ii,4})', stack_ACCR{ii,4}', max(stack_ACCR{ii,4})'];
        posAmps = [stack_ACCR{ii,1}', max((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,1})', stack_ACCR{ii,1}'];
        negAmps = [stack_ACCR{ii,1}', min((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,1})', stack_ACCR{ii,1}'];
        fill([posAmps,zeros(size(posAmps))+posAmps(1)],[Time,flip(Time)],'r','LineStyle','none')
        fill([negAmps,zeros(size(negAmps))+negAmps(1)],[Time,flip(Time)],'b','LineStyle','none')
        plot((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,4},'k')
    end
    for ii = 1:n_clust
        % GET TIMES FOR H-k-Vp solns
        Tpmp_pred = 2*sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2)*HSOLNS(ii,1);
        Tps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2 - rayp_pred.^2) - sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tppps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2) + sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tpsps_pred = 2*(sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2))*HSOLNS(ii,1);
        plot(rayp_pred,Tpmp_pred,'b','LineWidth',2)
        plot(rayp_pred,Tps_pred,'r','LineWidth',2)
        plot(rayp_pred,Tppps_pred,'r','LineWidth',2)
        plot(rayp_pred,Tpsps_pred,'b','LineWidth',2)
    end
else
    plot(rayp_pred,Tpmp_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tppps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tpsps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    for ii = 1:size(stack_ACCR,1)
        Time = [min(stack_ACCR{ii,4})', stack_ACCR{ii,4}', max(stack_ACCR{ii,4})'];
        posAmps = [stack_ACCR{ii,1}', max((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,1})', stack_ACCR{ii,1}'];
        negAmps = [stack_ACCR{ii,1}', min((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,1})', stack_ACCR{ii,1}'];
        fill([posAmps,zeros(size(posAmps))+posAmps(1)],[Time,flip(Time)],'r','LineStyle','none')
        fill([negAmps,zeros(size(negAmps))+negAmps(1)],[Time,flip(Time)],'b','LineStyle','none')
        plot((stack_ACCR{ii,COL}./max(abs(stack_ACCR{ii,COL}))).*scale+stack_ACCR{ii,1},stack_ACCR{ii,4},'k')
    end
    plot(rayp_pred,Tpmp_pred,'b','LineWidth',2)
    plot(rayp_pred,Tps_pred,'r','LineWidth',2)
    plot(rayp_pred,Tppps_pred,'r','LineWidth',2)
    plot(rayp_pred,Tpsps_pred,'b','LineWidth',2)
end
plot(rayp_pred,Tpmp_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tppps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tpsps_Hk,'k--','LineWidth',2)
axis([0.03 0.09 0 25])
axis ij
xlabel('Ray parameter (s/km)')


% RADIAL RECEIVER FUNCTIONS
subplot(3,3,3,'Position',[0.69,0.75,h_width,h_height])
bar(cell2mat(stack_ACCZ(:,1)),cell2mat(stack_ACCZ(:,3)))
max_y = max(cell2mat(stack_ACCZ(:,3)))+.2*max(cell2mat(stack_ACCZ(:,3)));
title('RADIAL RF')
set(gca,'xticklabel',{});
axis([0.03 0.09 0 max_y])
for ii = 1:size(stack_ACCZ,1)
    text(cell2mat(stack_ACCZ(ii,1)),cell2mat(stack_ACCZ(ii,3))+.1*max_y, ...
        sprintf('%0.0f',cell2mat(stack_ACCZ(ii,3))),'VerticalAlignment','top','HorizontalAlignment','center')
end

subplot(3,3,[6 9],'Position',[0.69,0.05,t_width,t_height])
scale = 0.013;
hold on
if n_clust > 1
    for ii = 1:n_clust
        plot(rayp_pred,Tpmp_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tppps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
        plot(rayp_pred,Tpsps_good(idx==ii,:),'Color',cc(ii,:),'LineWidth',0.5)
    end
    for ii = 1:size(stack_RFs,1)
        Time = [min(stack_RFs{ii,4})', stack_RFs{ii,4}', max(stack_RFs{ii,4})'];
        posAmps = [stack_RFs{ii,1}', max(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,1})', stack_RFs{ii,1}'];
        negAmps = [stack_RFs{ii,1}', min(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,1})', stack_RFs{ii,1}'];
        patch(posAmps,Time,[1 0 0],'EdgeColor', 'none')
        patch(negAmps,Time,[0 0 1],'EdgeColor', 'none')
        plot(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,4},'k')
    end
    for ii = 1:n_clust
        % GET TIMES FOR H-k-Vp solns
        Tpmp_pred = 2*sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2)*HSOLNS(ii,1);
        Tps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2 - rayp_pred.^2) - sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tppps_pred = (sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2) + sqrt(1/VPSOLNS(ii,1).^2 - rayp_pred.^2))*HSOLNS(ii,1);
        Tpsps_pred = 2*(sqrt(1/(VPSOLNS(ii,1)/VPVSSOLNS(ii,1)).^2-rayp_pred.^2))*HSOLNS(ii,1);
        plot(rayp_pred,Tpmp_pred,'Color',cc(ii,:),'LineWidth',2)
        plot(rayp_pred,Tps_pred,'Color',cc(ii,:),'LineWidth',2)
        plot(rayp_pred,Tppps_pred,'Color',cc(ii,:),'LineWidth',2)
        plot(rayp_pred,Tpsps_pred,'Color',cc(ii,:),'LineWidth',2)
    end
else
    plot(rayp_pred,Tpmp_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tppps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    plot(rayp_pred,Tpsps_good,'Color',[169/255 169/255 169/255],'LineWidth',0.5)
    for ii = 1:size(stack_RFs,1)
        Time = [min(stack_RFs{ii,4})', stack_RFs{ii,4}', max(stack_RFs{ii,4})'];
        posAmps = [stack_RFs{ii,1}', max(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,1})', stack_RFs{ii,1}'];
        negAmps = [stack_RFs{ii,1}', min(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,1})', stack_RFs{ii,1}'];
        patch(posAmps,Time,[1 0 0],'EdgeColor', 'none')
        patch(negAmps,Time,[0 0 1],'EdgeColor', 'none')
        plot(stack_RFs{ii,5}.*scale+stack_RFs{ii,1},stack_RFs{ii,4},'k')
    end
    plot(rayp_pred,Tpmp_pred,'r','LineWidth',2)
    plot(rayp_pred,Tps_pred,'r','LineWidth',2)
    plot(rayp_pred,Tppps_pred,'r','LineWidth',2)
    plot(rayp_pred,Tpsps_pred,'b','LineWidth',2)
end
plot(rayp_pred,Tpmp_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tppps_Hk,'k--','LineWidth',2)
plot(rayp_pred,Tpsps_Hk,'k--','LineWidth',2)
axis([0.03 0.09 0 25])
axis ij
box on
xlabel('Ray parameter (s/km)')

if plotme == 1
    print(char(fname),'-dpdf','-r0')
end

%% Find cluster solution



h=figure(4);
    % If H-k-Vp solution is good and errors are well quantified, click this
    GOOD = 0;
    PushButton1 = uicontrol(gcf,'Style','push','String','Good Solution','Position',[120 490 100 30]);
    set(PushButton1,'Callback',@PushB1_Callback);
    function PushB1_Callback(PushButton1,EventData)
        GOOD = get(PushButton1,'Value')
        if n_clust > 1
            clust_soln = cell2mat(inputdlg({'Solution Identification'},'Enter Solution ID',1,{'1'}));
            clust_soln = str2num(clust_soln);
            HkVp_RESULT(1) = H_good(Amps_good==max(Amps_good(idx==clust_soln)));
            HkVp_RESULT(2) = min(H_good(idx==clust_soln));
            HkVp_RESULT(3) = H_q(clust_soln,1);
            HkVp_RESULT(4) = H_q(clust_soln,2);
            HkVp_RESULT(5) = max(H_good(idx==clust_soln));
            HkVp_RESULT(6) = Vp_good(Amps_good==max(Amps_good(idx==clust_soln)));
            HkVp_RESULT(7) = min(Vp_good(idx==clust_soln));
            HkVp_RESULT(8) = Vp_q(clust_soln,1);
            HkVp_RESULT(9) = Vp_q(clust_soln,2);
            HkVp_RESULT(10) = max(Vp_good(idx==clust_soln));
            HkVp_RESULT(11) = VpVs_good(Amps_good==max(Amps_good(idx==clust_soln)));
            HkVp_RESULT(12) = min(VpVs_good(idx==clust_soln));
            HkVp_RESULT(13) = VpVs_q(clust_soln,1);
            HkVp_RESULT(14) = VpVs_q(clust_soln,2);
            HkVp_RESULT(15) = max(VpVs_good(idx==clust_soln));
            
            dlmwrite(fname2,[HkVp_RESULT(26), HkVp_RESULT(27),HkVp_RESULT(29),HkVp_RESULT(30),HkVp_RESULT(32),HkVp_RESULT(33),std_err_HkVp,std_err_Hk],'delimiter',' ')
            dlmwrite(fname2,[H_good Vp_good VpVs_good Amps_good idx],'delimiter',' ','-append');
            
%         else
%             HkVp_RESULT(1) = H_good(Amps_good==max(Amps_good));
%             HkVp_RESULT(2) = min(H_good);
%             HkVp_RESULT(3) = H_q(1,1);
%             HkVp_RESULT(4) = H_q(1,2);
%             HkVp_RESULT(5) = max(H_good);
%             HkVp_RESULT(6) = Vp_good(Amps_good==max(Amps_good));
%             HkVp_RESULT(7) = min(Vp_good);
%             HkVp_RESULT(8) = Vp_q(1,1);
%             HkVp_RESULT(9) = Vp_q(1,2);
%             HkVp_RESULT(10) = max(Vp_good);
%             HkVp_RESULT(11) = VpVs_good(Amps_good==max(Amps_good));
%             HkVp_RESULT(12) = min(VpVs_good);
%             HkVp_RESULT(13) = VpVs_q(1,1);
%             HkVp_RESULT(14) = VpVs_q(1,2);
%             HkVp_RESULT(15) = max(VpVs_good);            
        end
        close(3);close(4);close(5)
        %set(hpb1, 'position', [figstuff(3)/2 figstuff(4)/2 350 100]); %makes box bigger
    end

    % If solution or errors run into edges of gridsearch space, re-run
    % H-k-Vp with wider bounds
    newSP = [num2cell(minH) num2cell(maxH) num2cell(minVp) num2cell(maxVp) num2cell(minVpVs) num2cell(maxVpVs) num2cell(w(1)) num2cell(w(2)) num2cell(w(3))];
    PushButton2 = uicontrol(gcf,'Style','push','String','Re-run H-k-Vp stacking','Position',[300 490 120 30]);
    set(PushButton2,'Callback',@PushB2_Callback);
    function PushB2_Callback(PushButton2,EventData)
        %CLUST = get(PushButton2,'Value')
        prompt = {'min H:','max H:','min Vp:','max Vp:','min VpVs:','max VpVs:','w1 (Ps):','w2 (2P1S):','w3 (Pmp):'};
        dims = [1 20];
        definput = {num2str(minH),num2str(maxH),num2str(minVp),num2str(maxVp),num2str(minVpVs),num2str(maxVpVs) ...
            ,num2str(w(1)),num2str(w(2)),num2str(w(3))};
        newSP = inputdlg(prompt,'New H-k-Vp parameters',dims,definput);
        newSP = str2num(char(newSP))';
        
        HkVp_RESULT(26) = newSP(1);
        HkVp_RESULT(27) = newSP(2);
        HkVp_RESULT(29) = newSP(3);
        HkVp_RESULT(30) = newSP(4);
        HkVp_RESULT(32) = newSP(5);
        HkVp_RESULT(33) = newSP(6);
        HkVp_RESULT(36) = newSP(7)
        HkVp_RESULT(37) = newSP(8)
        HkVp_RESULT(38) = newSP(9)
        
        %msgbox(sprintf('%s \n%s','Cluster analysis needed.','Please close plot window.'));    
    end   

    % If multiple solutions exist in solution space, run cluster analysis
    PushButton3 = uicontrol(gcf,'Style','push','String','Run GMM Clustering','Position',[460 490 120 30]);
    set(PushButton3,'Callback',@PushB3_Callback);
    function PushB3_Callback(PushButton3,EventData)
        %CLUST = get(PushButton2,'Value')
        n_clust = cell2mat(inputdlg({'Enter cluster number'},'Number of Clusters',1,{'1'}));
        n_clust = str2num(n_clust);
        HkVp_RESULT(35) = n_clust
        %msgbox(sprintf('%s \n%s','Cluster analysis needed.','Please close plot window.'));    
    end    
       
        uiwait(h)
    


end

