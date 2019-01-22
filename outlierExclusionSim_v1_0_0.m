% Outlier identification simulation script
%
% Requires:         Palamedes (optional - for advanced adaptive methods)
%                   psignifit (optional - for psychometric curve fitting)
%
% Author(s):        Pete R Jones <petejonze@gmail.com>
% 
% Version History:  25/04/2016	PJ  Initial version
%                   27/04/2016	PJ  More methods
%                   29/04/2016	PJ  More methods
%                   13/05/2016	PJ  Cleaning & commenting
%                   31/05/2016	PJ  APP Submission
%                                               
%
% Copyright 2016 : P R Jones
% *********************************************************************
% 

    %%%%%%%%%
    %%  0  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        clc
        close all;
        clear all;
 

    %%%%%%%%%
    %%  1  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% User parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        tic();

        % check we have palamedes installed (if using Palamedes functions)
        % PAL_version();
        
        % general simulation params ---------------------------------------
        % N_DL_SIMS       = 10000
        % N_OUTLIER_SIMS  = 2000;
        % nTestPoints     = 32;
        % nOutliers       = 0:1:15;
        %
        % ultra quick version (debugging) [~1min]
        N_DL_SIMS       = 50;
        N_OUTLIER_SIMS  = 5;
        nTestPoints     = 32;
        nOutliers = [0 1 7 15];

        % simulated observer paramters ------------------------------------
        oparams             = [];
        oparams.pfunc       = 'logistic';
        oparams.mu          = {@truncnormrnd {[1 1], 8,3,8,30}};            % mean of posited gaussian CDF (psychometric function)  
        oparams.inoise      = {@truncnormrnd {[1 1], 2,2,2,15}};            % Determines observer sensitivity. Equal to slope parameter (standard deviation) when PsyFunc = Cumulative Gaussian
        oparams.lapseRate   = {@truncnormrnd {[1 1], 0.01,0.02,0,0.055}};   % [0 .05];  % chance of answering with a pure guess
        oparams.guessRate   = 0.5;                                          % chance of answering correctly by chance 

        % non-compliant (outlier) simulated observer paramters ------------
        oparams_poor             = oparams;
        oparams_poor.mu          = [15 20];
        oparams_poor.inoise      = [5 10];    % Determines observer sensitivity. Equal to slope parameter (standard deviation) when PsyFunc = Cumulative Gaussian
        oparams_poor.lapseRate   = [.5 .85];

        % set psychophysical paradigm parameters --------------------------
        pmethod         = 'mcs'; % 'QUEST, 'mcs'
        pparams         = [];
        pparams.minVal	= 1;
        pparams.maxVal 	= 64;
        cparams         = [];
        switch lower(pmethod)
            case 'quest'
                % QUEST ---------------------------------------------------
                pparams.priorAlphaRange = pparams.minVal:0.5:pparams.maxVal;
                pparams.beta            = 4;    % slope
                pparams.lambda          = 0.05; % estimate of lapse rate: the fraction of trials on which the observer presses blindly.
                pparams.gamma           = 0.05; % estimate of guess rate: the fraction of trials that will generate response 1 when intensity==-inf.
                pparams.PF              = @PAL_CumulativeNormal;
                pparams.stopCriterion   = 'reversals'; % 'trials';
                pparams.stopRule        = 8; % 64;
                % analysis
                cmethod                 = 'mean';
            case 'staircase'
                % Staircase -----------------------------------------------
                pparams.startVal       = 32;
                pparams.stepSize       = [4 2 1];
                pparams.downMod        = 1;
                pparams.nReversals     = [1 1 6];
                pparams.nUp            = [1 1 1];
                pparams.nDown          = [1 2 2];    % [50% 70.7% 70.7%]
                pparams.isAbsolute     = true; % false;
                pparams.minNTrials     = 0;
                pparams.maxNTrials     = inf;
                pparams.verbosity      = 0; % 2 to get something to watch
                % analysis
                cmethod                = 'reversals';
                cparams.nRevsToAvrg    = [2 4 -1]; % min 2, max 4, ignore 1st
            case 'mcs'
                % MCS -----------------------------------------------------
                pparams.stimLevels  = [1, 7, 13, 19, 25 31];
                pparams.nPerLevel   = 50;
                % analysis
                cmethod             = 'pfit';
                cparams.shape       = 'cumulative Gaussian'; % Gumbel; logistic; Weibull; linear
                cparams.n_intervals = 2;                     % 2AFC %correct plot
                cparams.nruns       = 100;                   % n bootstraps
                cparams.th_level    = .707;                  % .707 == the point targeted by 2down-1up staircase
                cparams.verbosity   = 0;
                cparams.plotMode    = 'no plot';
                cparams.lapseLim    = [];
            otherwise
                error('unknown method?');
        end
           

    %%%%%%%%%
    %%  2  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Simulate observer trial sequences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % simulate compliant-observers' trial sequences
        fprintf('\n\nSimulate trials given actual observer parameters (Compliant):\n');
        sims_compliant = PsychophysicalSimulator.generateTrialSequences(N_DL_SIMS, oparams, pmethod, pparams);
        DL_compliant = PsychophysicalSimulator.computeThresholds(sims_compliant, cmethod, cparams);

        % simulate non-compliant-observers' trial sequences (outliers)
        fprintf('\n\nSimulate trials given actual observer parameters (Non-Compliant):\n');
        sims_noncompliant = PsychophysicalSimulator.generateTrialSequences(N_DL_SIMS, oparams_poor, pmethod, pparams);
        DL_noncompliant = PsychophysicalSimulator.computeThresholds(sims_noncompliant, cmethod, cparams);
            
        % remove any NaN threshold estimates
        idx = isnan(DL_compliant);
        DL_compliant = DL_compliant(~idx);
        fprintf('\n\nRemoving %i NaN values from DL_compliant (%1.2f%%)\n', sum(idx), sum(idx)/length(idx));
        idx = isnan(DL_noncompliant);
        DL_noncompliant = DL_noncompliant(~idx);
        fprintf('Removing %i NaN values from DL_noncompliant (%1.2f%%)\n\n\n', sum(idx), sum(idx)/length(idx));
        
        
    %%%%%%%%%
    %%  3  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Check that compliant and outliers are actually distinguishable %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        figure('Position', [600 600 600 250]);
        set(0,'DefaultTextInterpreter','tex');

        % get data
        A = DL_compliant;
        B = DL_noncompliant;

        % Compute verlap coefficient
        nBins = 30;
        [~ , bins] = hist( [A; B] ,nBins);
        aCounts = hist( A , bins );
        bCounts = hist( B , bins );
        aP = aCounts./sum(aCounts);
        bP = bCounts./sum(bCounts);
        Bhattacharyya_Coefficient = sum(sqrt(aP.*bP));
        % alt measure - degree of overlap between histograms:
        minOfHists = min([aP; bP], [], 1);
        overlappedHist = sum(minOfHists);

        % plot
        hist(A,bins)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
        hold on;
        hist(B,bins)
        h1 = findobj(gca,'Type','patch');
        set(h1,'facealpha',0.75);

        % annotate
        txt = sprintf('BC = %1.2f', Bhattacharyya_Coefficient);
        title(txt);
        xlabel('DL');
        ylabel('N');

        % add legend
        legend(flipud(h1), 'Compliant', 'Non-compliant')

        
    %%%%%%%%%
    %%  4  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Apply heuristics to samples of actual data;                 %%%%%%%
    %%%  > run many simulations & systematically vary N outliers    %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % init output variables
        criterionTxt    = { 'SD(2)','SD(2.5)','SD(3)'                   ...
                            'rSD(2)','rSD(2.5)','rSD(3)'             	...
                            ,'Mix(0.05)','Mix(0.02)','Mix(0.01)'        ... 
                            ,'IQR(2)','IQR(2.5)','IQR(3)'               ... 
                         	,'prctile(95)','prctile(98)', 'prctile(99)' ...
                         	,'Tukey(1.5)', 'Tukey(3)'                   ...
                            ,'MAD$_{\textrm{n}}$(2)','MAD$_{\textrm{n}}$(2.5)','MAD$_{\textrm{n}}$(3)'         ...
                            ,'S$_{\textrm{n}}$(2)','S$_{\textrm{n}}$(2.5)','S$_{\textrm{n}}$(3)'};
        resp            = nan(nTestPoints, length(criterionTxt), N_OUTLIER_SIMS, length(nOutliers));
        targ            = nan(nTestPoints, length(criterionTxt), N_OUTLIER_SIMS, length(nOutliers));
        
        % run
        fprintf('*********** RUNNING OUTLIER SIMULATIONS ***********\n');    
        fprintf('Start date/time: %s\n---------------------------------------\n', datestr(now(),0))
        tic()
        warning off
        for l = 1:length(nOutliers)
            fprintf('\nN Outliers = %i of %i (%1.2f%%)\n', nOutliers(l), nTestPoints, nOutliers(l)/nTestPoints*100);
            for k = 1:N_OUTLIER_SIMS
                if mod(1:N_OUTLIER_SIMS, N_OUTLIER_SIMS/10)==0
                    fprintf('. ');
                end
                
                % get sample of data: mixture of compliant, and
                % non-compliant, using true (unknown) observer params
               	x = [randsample(DL_compliant,nTestPoints-nOutliers(l),false);
                    randsample(DL_noncompliant,nOutliers(l),false)];
                  
              	% 0. target, for scoring 
                targ(:,1,k,l) = [zeros(1,nTestPoints-nOutliers(l)) ones(1,nOutliers(l))];
                targ(:,:,k,l) = repmat(targ(:,1,k,l), [1 length(criterionTxt)]); % same for all

                % 1. SD ---------------------------------------------------
              	% apply "> mu + N SD" heuristics
                resp(:,1,k,l) = x > mean(x) + 2.0*std(x); % ALT: norminv(0.9772, mean(DL_test),std(DL_test))
                resp(:,2,k,l) = x > mean(x) + 2.5*std(x);
                resp(:,3,k,l) = x > mean(x) + 3.0*std(x);
                
                % 2. recursive SD -----------------------------------------
                maxN = 3;
                lambda = [2 2.5 3];
                for j = 1:3
                    xx = x;
                    for i = 1:maxN
                        idx = xx>(nanmean(xx)+nanstd(xx)*lambda(j));
                        if sum(idx)==0
                            break
                        end
                        xx(idx) = NaN;
                    end
                    resp(:,3+j,k,l) = isnan(xx);
                end
                
                % 3. mixture model ----------------------------------------
                % define mixture shapes
                fC_name = 'gamma';
                fO_name = 'normal';
                % define mixture function form
                fS = @(x, w, fC_param1,fC_param2, fO_param1,fO_param2) w*pdf(ProbDistUnivParam(fC_name, [fC_param1 fC_param2]), x) + (1-w)*pdf(ProbDistUnivParam(fO_name, [fO_param1 fO_param2]), x);
              	% fit non-mixture model (just compliant, fC, component)
                simpleFit = fitdist(x, fC_name);
                % starting guesses
                w_start  	= .95;
                fC_start	= simpleFit.Params;
                fO_start	= [max(x) simpleFit.Params(2)/2];
                start       = [w_start fC_start fO_start];
                % parameter bounds
                lb = [.66 -Inf 0.01 min(simpleFit.icdf(.75),max(x))   simpleFit.Params(2)/10];
                ub = [1    Inf  Inf Inf                                     Inf];
                % n iterations
                options = statset('MaxIter',999, 'MaxFunEvals',999);
                % run! find MLE mixture model
                paramEsts = mle(x, 'pdf',fS, 'start',start, ...
                   'lower',lb, 'upper',ub, 'options',options);
                [w_est,fC_est,fO_est] = deal(paramEsts(1),paramEsts(2:3),paramEsts(4:5));
                % establish pdfs
                fC = ProbDistUnivParam(fC_name, fC_est);
                fO = ProbDistUnivParam(fO_name, fO_est);
                % post-hoc exclusion of the mixture component (i.e., if
                % compliant/outlier distributions too overlapping)
                % compute overlap
                xFit = linspace(min(x),max(x),1000);
                d1 =     w_est*fC.pdf(xFit);
                d2 = (1-w_est)*fO.pdf(xFit);
                [~, index] = min(abs(d1./d2-1)); % http://stackoverflow.com/questions/29509701/how-to-find-intersection-of-two-distribution-in-matlab
                % determine exclusion
                if fC.cdf(xFit(index))<.99
                   fC = simpleFit;
                end
                % establish exclusion criterion based on just the compliant
                % component of the mixture model
                resp(:,7,k,l) = x > fC.icdf(1-0.0228); % alpha = 0.05, equiv to 2SD
                resp(:,8,k,l) = x > fC.icdf(1-0.01); % alpha = 0.01
                resp(:,9,k,l) = x > fC.icdf(1-0.001); % alpha = 0.001, equiv to 3SD

                % 4. IQR --------------------------------------------------
                % apply "> median + N iqr" heuristics
                resp(:,10,k,l) = x > median(x) + 2*iqr(x);
                resp(:,11,k,l) = x > median(x) + 2.5*iqr(x);
                resp(:,12,k,l) = x > median(x) + 3*iqr(x);

                % 5. prctile ----------------------------------------------
                resp(:,13,k,l) = x > prctile(x, 95);
                resp(:,14,k,l) = x > prctile(x, 98);
                resp(:,15,k,l) = x > prctile(x, 99);

                % 6. Tukey ------------------------------------------------
                % apply Tukey's "> Q3 + N iqr" [boxplot/fence] heuristics
                resp(:,16,k,l) = x > prctile(x, 75) + 1.5*iqr(x);
                resp(:,17,k,l) = x > prctile(x, 75) + 3*iqr(x);
                
                % 7. MAD_n ------------------------------------------------
                MADn = 1.4826*mad(x,1);
                resp(:,18,k,l) = (x-median(x))/MADn > 2.0;
                resp(:,19,k,l) = (x-median(x))/MADn > 2.5;
                resp(:,20,k,l) = (x-median(x))/MADn > 3.0;
                
                % 8. S_n --------------------------------------------------
                [Sn, x_j] = RousseeuwCrouxSn(x);
                resp(:,21,k,l) = x_j/Sn > 2.0;
                resp(:,22,k,l) = x_j/Sn > 2.5;
                resp(:,23,k,l) = x_j/Sn > 3.0;
            end
        end
        fprintf('\n done!\n');
        warning off
        toc();
        fprintf('-----------------------------------------------\n');

        
    %%%%%%%%%
    %%  5  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Score each method                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        % score individual simulations
        nHit      	= permute(sum(resp==1 & targ==1), [3 2 4 1]);
        nFalseAlarm	= permute(sum(resp==1 & targ==0), [3 2 4 1]);
        nMiss     	= permute(sum(resp==0 & targ==1), [3 2 4 1]);
        nCorrectRej	= permute(sum(resp==0 & targ==0), [3 2 4 1]);
        pCorrect   	= permute(mean(resp==targ),       [3 2 4 1]);
 
    	% tabulate means of all simulation, for plotting

        hit_N = reshape(mean(nHit),length(criterionTxt),length(nOutliers))';
        hit_Rate = bsxfun(@rdivide,hit_N',nOutliers)';
        FA_N = reshape(mean(nFalseAlarm),length(criterionTxt),length(nOutliers))';
        FA_Rate = bsxfun(@rdivide,FA_N',nTestPoints-nOutliers)';
        pC = reshape(mean(pCorrect),length(criterionTxt),length(nOutliers))';
        
    	% tabulate StD of all simulation, for plotting
        hit_N_SD = reshape(std(nHit),length(criterionTxt),length(nOutliers))';
        FA_N_SD = reshape(std(nFalseAlarm),length(criterionTxt),length(nOutliers))';
        pC_SD = reshape(std(pCorrect),length(criterionTxt),length(nOutliers))';
        
        
    %%%%%%%%%
    %%  6  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save data                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       fprintf('saving..\n')
       fn = sprintf('outlierExclusionSim_v0_1_2__nDL=%i_NSIMS=%i_nTest=%i_NnOut=%i__%s.mat',N_DL_SIMS,N_OUTLIER_SIMS,nTestPoints,length(nOutliers),datestr(now(),30));
       save(fn);
       fprintf('  done!\n')

       
    %%%%%%%%%
    %%  7  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot results                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        %      -Std-  -rSD-  -Mix-  -IQR-  -Prc-  -T-  -MAD-  -Sn--
        idx = [1 0 1  0 1 0  0 0 1  1 0 0  1 0 0  1 0  0 0 1  0 0 1]==1;
        
        % Plot MEAN performance -------------------------------------------
 
        figure() % Rates
        %1. Hits
        subplot(3,1,1);
        hDat = plot(nOutliers, hit_Rate(:,idx), '-');
        ylabel('Hit Rate');
        set(hDat(end-1:end),'linestyle',':');
        %2. False Alarms
        subplot(3,1,2);
        hDat = plot(nOutliers, FA_Rate(:,idx),'-');
        ylabel('FA Rate');
        set(hDat(end-1:end),'linestyle',':');
        %3. Percent Correct
        subplot(3,1,3);
        hDat = plot(nOutliers, pC(:,idx),'-');
        ylabel('\%Correct');
        set(hDat(end-1:end),'linestyle',':');
       	% add legend
        hLeg = legend(criterionTxt(idx));
        set(hLeg, 'Location','East');