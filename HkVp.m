function [ tern_info, HkVp_RESULT, HkVp_good, Hk_good, std_err_HkVp, std_err_Hk ] = HkVp(stack_ACCZ,stack_ACCR,stack_RFs,SP,rayp,zonly,AssVp,wtog,ptog,w,COL,fname2,fname21,fname3)
% This function calculates the travel times for phases in receiver
% functions and autocorrelations to be used later to do H-k-Vp stacking
  % stack_* are binned and stacked waveforms from EQ_AUTOCORR    
  % H, Vp, VpVs are vectors that contain the gridsearch range for that
  %   variable. 
  % rayp is a vector of the center of bins from stacking (obtained in
  %   previous script EQ_AUTOCORR)
  % zonly specifies whether you only want to use vertical ACs (1) or both
  %   vertical and radial ACs
  % AssVp specifies assumed Vp for H-k stacking
  % wtog is toggle that specifies if you want to weight by bin population (1=yes, 0=no)
  % ptog is toggle that specifies if you want to discount amplitudes for
  %   predicted arrivals if they have the wrong polarity (1=yes, 0=no)
  % w is weighting vector of different phases
  % COL define whether linear or pw stacks are plotted (5 = linear, 6 = pwstack)
  % fname2 is name of matfile containing "good" solutions
  % fname21 is name of textfile containing "good" solutions
  % fname3 is textfile with results of individual station


%% Get timing of arrivals

H = SP(1):SP(3):SP(2);
Vp = SP(4):SP(6):SP(5);
VpVs = SP(7):SP(9):SP(8);

c=0;
Tpmp = zeros(length(H)*length(Vp)*length(VpVs),length(rayp));
Tps = zeros(length(H)*length(Vp)*length(VpVs),length(rayp));
Tpps = zeros(length(H)*length(Vp)*length(VpVs),length(rayp));
H_vect = zeros(length(H)*length(Vp)*length(VpVs),1);
Vp_vect = zeros(length(H)*length(Vp)*length(VpVs),1);
VpVs_vect = zeros(length(H)*length(Vp)*length(VpVs),1);

for j = 1:length(H) 
    for k = 1:length(Vp)
        for l = 1:length(VpVs)
            c=c+1;
            H_vect(c) = H(j);
            Vp_vect(c) = Vp(k);
            VpVs_vect(c) = VpVs(l);
            Vs = Vp(k)./VpVs(l);
            Tpmp(c,:) = 2.*H(j).*sqrt(1/Vp(k).^2 - rayp.^2); % Tpmp
            Tps(c,:) = H(j).*(sqrt(1/Vs.^2 - rayp.^2) - sqrt(1/Vp(k).^2-rayp.^2)); % Tps
            Tpps(c,:) = H(j).*(sqrt(1/Vs.^2 - rayp.^2) + sqrt(1/Vp(k).^2-rayp.^2));  % Tpps 
        end
    end
end



% Extract populated bin stacks to find relevant arrivals
rayp2 = cell2mat(stack_ACCZ(:,1));
rayp_pop = ismember(rayp(:),round(rayp2'*10^4)/10^4)'; % check which bins are populated
newTps = Tps(:,rayp_pop);
newTpmp = Tpmp(:,rayp_pop);
newTpps = Tpps(:,rayp_pop);
weightsACC = cell2mat(stack_ACCZ(:,3));
nevents = sum(weightsACC);

%% Perform H-k-Vp and H-k Stacking
for m = 1:length(rayp2)
    % Get amplitudes of Pmp phase from autocorrelations
    % If wtog == 1, contributions of individual bins scaled to number
    % of events in that bin. If wtog == 0, all bins are given equal
    % weighting to the final stack
    if wtog == 1
        %   Note: weightsACC is same for RFs if you use the traces that made
        %   the RFs for autocorrelation
        AmpZ(:,m) = interp1(stack_ACCZ{m,4},stack_ACCZ{m,COL},newTpmp(:,m))*weightsACC(m);
        AmpR(:,m) = interp1(stack_ACCR{m,4},stack_ACCR{m,COL},newTpmp(:,m))*weightsACC(m);
        Amp_Ps(:,m) = interp1(stack_RFs{m,4},stack_RFs{m,5},newTps(:,m)*weightsACC(m));
        Amp_Pps(:,m) = interp1(stack_RFs{m,4},stack_RFs{m,5},newTpps(:,m)*weightsACC(m)); 
    else
        AmpZ(:,m) = interp1(stack_ACCZ{m,4},stack_ACCZ{m,COL},newTpmp(:,m));
        AmpR(:,m) = interp1(stack_ACCR{m,4},stack_ACCR{m,COL},newTpmp(:,m));    
        Amp_Ps(:,m) = interp1(stack_RFs{m,4},stack_RFs{m,5},newTps(:,m));
        Amp_Pps(:,m) = interp1(stack_RFs{m,4},stack_RFs{m,5},newTpps(:,m)); 
    end
end        


if wtog == 1
    if ptog == 1
        % If you only want to consider arrivals with the expected
        % polarity of the arrival (amplitudes that have the opposite
        % polarity of the expected phase will not contribute to
        % stacking)
        AmpZ(AmpZ>0) = 0;
        AmpR(AmpR>0) = 0;
        Amp_Ps(Amp_Ps<0) = 0;
        Amp_Pps(Amp_Pps<0) = 0;
        sumAmp_Zpmp = sum(AmpZ./sum(weightsACC),2);
        sumAmp_Rpmp = sum(AmpR./sum(weightsACC),2);
        sumAmp_Ps = sum(Amp_Ps./sum(weightsACC),2);
        sumAmp_Pps = sum(Amp_Pps./sum(weightsACC),2); 
    else
        % Sum all amplitudes, even if polarity is opposite what is
        % expected (could lead to destructive stacking)
        sumAmp_Zpmp = sum(AmpZ./sum(weightsACC),2);
        sumAmp_Rpmp = sum(AmpR./sum(weightsACC),2);
        sumAmp_Ps = sum(Amp_Ps./sum(weightsACC),2);
        sumAmp_Pps = sum(Amp_Pps./sum(weightsACC),2); 
    end
else
    % No normalization by events in bin. All bins contribute equally to
    % stacked solution
    if ptog == 1
        % For ptog description, look up ~20 lines
        AmpZ(AmpZ>0) = 0;
        AmpR(AmpR>0) = 0;
        Amp_Ps(Amp_Ps<0) = 0;
        Amp_Pps(Amp_Pps<0) = 0;
        sumAmp_Zpmp = sum(AmpZ,2);
        sumAmp_Rpmp = sum(AmpR,2);
        sumAmp_Ps = sum(Amp_Ps,2);
        sumAmp_Pps = sum(Amp_Pps,2); 
    else
        sumAmp_Zpmp = sum(AmpZ,2);
        sumAmp_Rpmp = sum(AmpR,2);
        sumAmp_Ps = sum(Amp_Ps,2);
        sumAmp_Pps = sum(Amp_Pps,2); 
    end
end

  
% Make matrix of all results
result_mat = [ H_vect Vp_vect VpVs_vect sumAmp_Zpmp sumAmp_Rpmp sumAmp_Ps sumAmp_Pps ];  

% Calculate ACC weight using either 0 (sum Z and R ACC) or 1 (only use
% Z ACC), and normalize to largest (negative) arrival but do not change
% sign
if zonly == 2
    weight_ACC = (result_mat(:,4)/abs(min(result_mat(:,4))) + result_mat(:,5)/abs(min(result_mat(:,5))))./2;
elseif zonly == 1
    weight_ACC = result_mat(:,4)/abs(min(result_mat(:,4)));
end

% Normalize Ps and Pps arrivals
weight_Ps = result_mat(:,6)./max(result_mat(:,6));
weight_Pps =result_mat(:,7)./max(result_mat(:,7));

% Calculate stacked amplitudes
sum_weight_HkVp = -w(3).*weight_ACC + w(1).*weight_Ps + w(2).*weight_Pps;



%% Find H-Vp-k solutions
    % tern_info contains H-k-Vp gridsearch results
    tern_info = [H_vect Vp_vect VpVs_vect sum_weight_HkVp];
    tern_info(:,4) = tern_info(:,4)/max(tern_info(:,4)); % Normalize so max amplitude = 1 
    
    % Calculated standard error and subtract from maximum amplitude in stack
    % to constrain data in "good" solutions (Eaton et al., 2006)
    std_err_HkVp = max(tern_info(:,4)) - sqrt(std(tern_info(:,4)).^2./nevents);

    % Find all results withing standard error of max in HkVp stack
    H_good = tern_info(tern_info(:,4)>=std_err_HkVp,1);
    Vp_good = tern_info(tern_info(:,4)>=std_err_HkVp,2);
    VpVs_good = tern_info(tern_info(:,4)>=std_err_HkVp,3);
    Amps_good = tern_info(tern_info(:,4)>=std_err_HkVp,4);
    
    HkVp_good = [H_good Vp_good VpVs_good Amps_good];
    % Write out file of "good" solutions (for GMT plotting)
    %dlmwrite(fname2,[min(H), max(H),min(Vp), max(Vp), min(VpVs), max(VpVs) std_err_HkVp],'delimiter',' ')
    %dlmwrite(fname2,[H_good Vp_good,VpVs_good,Amps_good],'delimiter',' ','-append')    


    %% Find H-K solutions
    indhk = find(round(Vp_vect*100)/100==AssVp);
    Hkmat = result_mat(indhk,[1 3 6 7]);
    
    % Normalize amplitudes for each phase, define weights based on first
    % two indices of weighting vector, create H-k stack
    weight_Ps = Hkmat(:,3)./max(Hkmat(:,3));
    weight_Pps =Hkmat(:,4)./max(Hkmat(:,4));
    sum_weight_Hk = (w(1)/(w(1)+w(2))).*weight_Ps + (w(2)/(w(1)+w(2))).*weight_Pps;
    sum_weight_Hk = sum_weight_Hk./max(sum_weight_Hk);
    
    % Find standard error of Hk stack
    std_err_Hk = max(sum_weight_Hk) - sqrt(std(sum_weight_Hk).^2./nevents);

    % Get largest amplitude from Hk solution and "good" solutions
    H_Hkbest = Hkmat(sum_weight_Hk==1,1);
    VpVs_Hkbest = Hkmat(sum_weight_Hk==1,2);
    Amps_Hkgood = sum_weight_Hk(sum_weight_Hk>=std_err_Hk);
    H_Hkgood = Hkmat(sum_weight_Hk>=std_err_Hk,1);
    VpVs_Hkgood = Hkmat(sum_weight_Hk>=std_err_Hk,2);
    
    Hk_good = [H_Hkgood VpVs_Hkgood Amps_Hkgood];

    %% Calculate Errors: interquartile distribution and min/max of value for box and whisker plots

    % Get best solution
    H_best = tern_info(tern_info(:,4)==1,1);
    Vp_best = tern_info(tern_info(:,4)==1,2);
    VpVs_best = tern_info(tern_info(:,4)==1,3);

    % Get errors for HkVp and Hk based on population of "good" solutions
    H_q = quantile(H_good,[.159 .841]); % Find top 68.2% of H (equiv to 1 standard dev for normally distributed data)
    Vp_q = quantile(Vp_good,[0.159 0.841]);
    VpVs_q = quantile(VpVs_good,[0.159 0.841]);
    
    H_Hkq = quantile(H_Hkgood,[0.159 0.841]);
    VpVs_Hkq = quantile(VpVs_Hkgood,[0.159 0.841]);
    %old wrong number used [.196 .814]
    
    %% Array-itize the results
    %HkVp_RESULT = [H_best H_q(1) H_q(2) Vp_best Vp_q(1) Vp_q(2) VpVs_best VpVs_q(1) VpVs_q(2) H_Hkbest H_Hkq(1) H_Hkq(2) VpVs_Hkbest VpVs_Hkq(1) VpVs_Hkq(2) ...
        %min(H) max(H) step_H min(Vp) max(Vp) step_Vp min(VpVs) max(VpVs) step_VpVs];
    HkVp_RESULT = [H_best min(H_good) H_q(1) H_q(2) max(H_good) Vp_best min(Vp_good) Vp_q(1) Vp_q(2) max(Vp_good) ...
        VpVs_best min(VpVs_good) VpVs_q(1) VpVs_q(2) max(VpVs_good) H_Hkbest min(H_Hkgood) H_Hkq(1) H_Hkq(2) max(H_Hkgood) ...
        VpVs_Hkbest min(VpVs_Hkgood) VpVs_Hkq(1) VpVs_Hkq(2) max(VpVs_Hkgood), SP(1), SP(2), SP(3), SP(4), SP(5), SP(6), SP(7), SP(8), SP(9) 1 ...
        w(1) w(2) w(3)];
    % Save out *_terninfo.mat with tern_info data and pdf parameters
    %BHinfo = [H_best min(H_good) H_q(1) H_q(2) max(H_good) Vp_best min(Vp_good) Vp_q(1) Vp_q(2) max(Vp_good) ...
        %VpVs_best min(VpVs_good) VpVs_q(1) VpVs_q(2) max(VpVs_good) H_Hkbest min(H_Hkgood) H_Hkq(1) H_Hkq(2) max(H_Hkgood) ...
        %VpVs_Hkbest min(VpVs_Hkgood) VpVs_Hkq(1) VpVs_Hkq(2) max(VpVs_Hkgood)];
    
    % Write out file of "good" solutions (for GMT plotting)
    dlmwrite(fname2,[min(H), max(H),min(Vp), max(Vp), min(VpVs), max(VpVs) std_err_HkVp],'delimiter',' ')
    dlmwrite(fname2,[H_good Vp_good,VpVs_good,Amps_good],'delimiter',' ','-append')
    save(fname21,'tern_info','HkVp_RESULT','std_err_HkVp','std_err_Hk','HkVp_good','Hk_good','AssVp');
    % write out file with box and whisker plot format for GMT
    dlmwrite(fname3,[H_best min(H_good) H_q(1) H_q(2) max(H_good) Vp_best min(Vp_good) Vp_q(1) Vp_q(2) max(Vp_good) ...
        VpVs_best min(VpVs_good) VpVs_q(1) VpVs_q(2) max(VpVs_good) H_Hkbest min(H_Hkgood) H_Hkq(1) H_Hkq(2) max(H_Hkgood) ...
        VpVs_Hkbest min(VpVs_Hkgood) VpVs_Hkq(1) VpVs_Hkq(2) max(VpVs_Hkgood)],'delimiter',' ');
    
end

