function status = checkingProgression(t, ~, flag)
    %CHECKINGPROGRESSION Display progress during ODE integration
    %  This function is used as an OutputFcn for ode45 to show a progress
    %  bar during the numerical integration.
    %
    %  Inputs:
    %     t    - Current time (scalar or vector depending on flag)
    %     y    - Current state (not used but required by MATLAB interface)
    %     flag - String indicating the call context: 'init', [], or 'done'
    %
    %  Output:
    %     status - 0 to continue integration, 1 to stop
    
    persistent count t_start t_end;
    
    if strcmp(flag, 'init')
        % Initialization phase: called once at the beginning
        count = 0;
        t_start = t(1);
        t_end = t(end);  % t contains [t_start, t_end] during init
        fprintf('Starting integration from t = %.2f to t = %.2f\n', t_start, t_end);
        
    elseif isempty(flag)
        % Integration phase: called at each successful step
        count = count + 1;
        
        % Display progress every 100 evaluations
        if mod(count, 100) == 0
            t_current = t(end);  % if t is a vector, take last value
            progress = (t_current - t_start) / (t_end - t_start);
            percent = progress * 100;
            
            % Create visual progress bar (40 characters wide)
            num_bars = round(progress * 40);  % Number of '=' characters
            bar = ['[' repmat('▄', 1, num_bars) repmat(' ', 1, 40-num_bars) ']'];
            
            % Update on same line using \r (carriage return)
            fprintf('\r%s %5.1f%% | t = %8.4f | iter: %d', bar, percent, t_current, count);
        end
        
    elseif strcmp(flag, 'done')
        % Finalization phase: called once at the end
        fprintf('\n'); % New line to preserve final progress bar
        fprintf('Integration completed! Total evaluations: %d\n', count);
        
        % Reset persistent variables for next call
        count = [];
        t_start = [];
        t_end = [];
    end
    
    status = 0;  % Always continue 
end