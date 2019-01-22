classdef PsychophysicalSimulator < handle
%#####
%
%   #####
%
% @Requires:        Palamedes toolkit
%
% @Example:         ######
%
% @See also:        AdaptiveTrack_example.m
%
% @Earliest compatible Matlab version:	v2008
%
% @Author:          Pete R Jones
%
% @Current Verion:  1.0.0
% @Version History: v1.0.0	PJ 27/03/2016    Initial build.
    
    
    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================
    
    methods(Static)
        
        function sims = generateTrialSequences(nSims, oparams, pmethod, pparams)
            fprintf('-----------------------------------------------\n');
            fprintf('Simulating observer trial sequences');
            tic();
            
            % init output variables -------------------------------------------
            sims = [];
            switch lower(pmethod)
                case 'staircase'        
                    sims.delta      = cell(nSims, 1);
                    sims.anscorrect	= cell(nSims, 1);
                    sims.reversals	= cell(nSims, 1);
                case 'quest'
                    sims.delta      = cell(nSims, 1);
                    sims.anscorrect	= cell(nSims, 1);
                    sims.mean       = nan(nSims, 1);
                    sims.mode       = nan(nSims, 1);
                case 'mcs'
                    sims.delta      = cell(nSims, 1);
                    sims.anscorrect	= cell(nSims, 1);
                otherwise
                    error('Unknown psychophysical paradigm: %s', pmethod);
            end
        
                        
            % RUN -------------------------------------------------------------
            for i = 1:nSims
                
                % print progress update(s) to the monitor
                if ismember(i, round(linspace(1,nSims,11)))
                    fprintf('. ')
                end
                
                % give estimate of time remaining
                if i==2 % ignore first run, as will be caching functions
                    startTime = tic();
                elseif i==12
                    cycleTime = toc(startTime)/10;
                    fprintf('\n  > cycle time is %1.2f secs per Sim (approx time remaining: %1.2f mins)\n  ', cycleTime, ((nSims-11)*cycleTime)/60)
                end
                
                % generate a new observer
                O = Observer(oparams);
                
                % initialise psychophysical paradigm & generate trial
                % sequences and use them to estimate threshold
                switch lower(pmethod)
                    case 'staircase'
                        % initialise
                        P = AdaptiveTrack(pparams);
                        % run
                        while ~P.isFinished()
                            % Get recommended stimulus level
                            x = P.getDelta();
                            % Simulate response
                            anscorrect = O.getResponse(x);
                            % Update algorithm
                            P.update(anscorrect);
                        end
                        % store raw data
                        sims.delta{i}     	= P.deltaHistory();
                        sims.anscorrect{i}  = P.wasCorrectHistory();
                        sims.reversals{i}	= P.getReversals();
                    case 'quest'
                        P = PAL_AMRF_setupRF(PAL_AMRF_setupRF(), 'priorAlphaRange',pparams.priorAlphaRange, 'beta',pparams.beta, 'lambda',pparams.lambda, 'gamma',pparams.gamma, 'PF',pparams.PF, 'stopRule',pparams.stopRule);
                        while ~P.stop()
                            x = P.xCurrent();
                            anscorrect = O.getResponse(x);
                            P = PAL_AMRF_updateRF(P, x, anscorrect);
                        end
                        % store raw data
                        sims.delta{i}     	= P.x;
                        sims.anscorrect{i}  = P.response;
                        sims.mean(i)        = P.mean;
                        sims.mode(i)        = P.mode;
                    case 'mcs'
                        % initialise
                        P = ConstantStimuli(pparams.stimLevels, pparams.nPerLevel);
                        % run
                        while ~P.isFinished()
                            % Get recommended stimulus level
                            x = P.getDelta();
                            % Simulate response
                            anscorrect = O.getResponse(x);
                            % Update algorithm
                            P.update(anscorrect);
                        end
                        % store raw data
                        sims.delta{i}     	= P.trialval();
                        sims.anscorrect{i}  = P.anscorrect();
                    otherwise
                        error('Unknown psychophysical paradigm: %s', pmethod);
                end
                sims.pmethod = lower(pmethod);
            end % end of N Sims
           
            fprintf(' done!\n');
            toc();
            fprintf('-----------------------------------------------\n');
        end
        
        function DL = computeThresholds(sims, cmethod, cparams)
            fprintf('-----------------------------------------------\n');   
            fprintf('Computing DLs..\n');
            tic();
            
            % init output variables -------------------------------------------
            nSims   = length(sims.delta);
            DL      = nan(nSims, 1);
            
            % validate
            switch cmethod
                case 'reversals'
                    if ~any(strcmpi(sims.pmethod, 'staircase'))
                        error('only makes sense to compute threshold by averaging reversals when using a staircase');
                    end
                case 'mean'
                    if ~any(strcmpi(sims.pmethod, 'quest'))
                        error('only makes sense to compute threshold from posterior MEAN when using quest');
                    end
                case 'mode'
                    if ~any(strcmpi(sims.pmethod, 'quest'))
                        error('only makes sense to compute threshold from posterior MODE when using quest');
                    end
                case 'pfit'
                    % no checks
                otherwise
                    error('Unknown DL compution paradigm: %s', cmethod);
            end
                    
                        
            % RUN -------------------------------------------------------------
            for i = 1:nSims
                switch cmethod
                    case 'reversals'
                        % get data
                        vals        = sims.reversals{i};
                        % compute n reversals available
                        nvals = length(vals);
                        if length(cparams.nRevsToAvrg)>2
                            nvals = nvals + cparams.nRevsToAvrg(3); % subtract first N reversals
                        end
                       	maxNEvenReversals = nvals - (mod(nvals,2)==1); % get max even number of (remaining) reversals
                        % constrain to specified range
                        if maxNEvenReversals < cparams.nRevsToAvrg(1)
                            warning('insufficient reversals (%i) to reach specified lower bound (%i). DL set to NaN', maxNEvenReversals, cparams.nRevsToAvrg(1));
                        else
                            nRevsToAvrg = min(maxNEvenReversals, cparams.nRevsToAvrg(2));
                            % compute
                            DL(i)       = mean(vals(end-(nRevsToAvrg-1):end));
                        end
                    case 'mean'
                        DL(i)       = sims.mean(i);
                    case 'mode'
                        DL(i)       = sims.mode(i);
                    case 'pfit'
                        % fit psychometric function
                        xyn = [sims.delta{i} sims.anscorrect{i} ones(size(sims.anscorrect{i}))];
                        fit = pfit(xyn, 'n_intervals',cparams.n_intervals, 'verbose',cparams.verbosity, 'runs',cparams.nruns, 'shape',cparams.shape, 'plot_opt',cparams.plotMode, 'lambda_limits',cparams.lapseLim);
                        % get threshold
                        DL(i) = findthreshold(cparams.shape, fit.params.est, cparams.th_level, 'performance');
                    otherwise
                        error('Unknown DL compution paradigm: %s', cmethod);
                end
                              
            end % end of N Sims
            
            fprintf(' done!\n');
            toc();
            fprintf('-----------------------------------------------\n');
        end
    end
    
end