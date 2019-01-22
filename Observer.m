classdef Observer < handle
    %#####
    %
    %   ####
    %
    % @Requires the following toolkits: <none>
    %
    % @Constructor Parameters:
    %
    %     	######
    %
    %
    % @Example:         O = Observer(Observer.getDummyParams())
    %                   anscorrect = Observer.getResponse(stimLevel);
    %
    % @See also:        #####
    %
    % @Earliest compatible Matlab version:	v2008
    %
    % @Author:          Pete R Jones
    %
    % @Creation Date:	17/03/16
    % @Last Update:     17/03/16
    %
    % @Current Verion:  1.0.0
    % @Version History: v1.0.0	PJ 17/03/16    Initial build.
    %
    % @Todo:            lots
    
    properties (GetAccess = 'public', SetAccess = 'private')
        % user specified parameters
        pfunc;
        mu;
        inoise;
        lapseRate;
        guessRate;
    end
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
    
    methods (Access = 'public')
        
        %% == CONSTRUCTOR =================================================
        
        function obj=Observer(params)
            % parse inputs & set specified parameter values ---------------
            obj.pfunc     	= params.pfunc;
            %if ~any(strcmpi(params.pfunc, {'normal', 'logistic'}))
            %    error('not yet supported');
            %end
            obj.mu          = obj.parseInputParam(params.mu);
            obj.inoise      = obj.parseInputParam(params.inoise);
            obj.lapseRate   = obj.parseInputParam(params.lapseRate);
            obj.guessRate   = obj.parseInputParam(params.guessRate);
            % try running to test
            cdf(obj.pfunc, 1, obj.mu, obj.inoise);
        end
        % Destructor
        function obj = delete(obj)        
            clear obj;
        end
        
        %% == METHODS =====================================================
        
        function anscorrect = getResponse(obj, stimLevel)
            
            % Simulate response
            if rand() < obj.lapseRate
                anscorrect = rand()> obj.guessRate;
            else
                %anscorrect = stimLevel > normrnd(obj.mu, obj.inoise);
                anscorrect = stimLevel > random(obj.pfunc, obj.mu, obj.inoise);
            end
            
        end
        
        
        
    end
    
    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================
    
    methods(Static)
        % useful when debugging
        function params = getDummyParams()
            params = [];
            params.pfunc = 'normal';
            params.mu = [10 15];
            params.inoise = [3 6];
            params.lapseRate = [0 .05]; % (truncnormrnd([1 1], 0.01, 0.0125, 0, 0.05)
            params.guessRate = .5;
        end
    end
    
    %% ====================================================================
    %  -----PRIVATE METHODS-----
    %$ ====================================================================
    
    methods(Access = 'private')

        function paramOut = parseInputParam(obj, paramIn)
            if iscell(paramIn)
                paramOut = paramIn{1}(paramIn{2}{:}); %#ok<*PROP> % e.g., paramIn = {@truncnormrnd {[1 1] 10 3 8 12}}
            elseif ismatrix(paramIn)
                if numel(paramIn)==1
                    paramOut = paramIn;
                elseif numel(paramIn)==2
                    paramOut = unifrnd(paramIn(1), paramIn(2));
                elseif numel(paramIn)==4
                    paramOut = truncnormrnd([1 1], paramIn(1), paramIn(2), paramIn(3), paramIn(4));
                else
                    error('unknown input length');
                end
            else
                error('unknown input type');
            end
        end
    end
end