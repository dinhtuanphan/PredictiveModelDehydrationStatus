function [plateau_start_index, plateau_end_index, plateau_true]...
    = A_plateau(time, data_value, time_start, time_end, time_window_min, threshold, direction)

    % time: time data, data_value: data
    % time_start: the start time of analysis, time_end: the end time of analysis
    % threshold: mM / min typical value = 0.1 mM/min
    % direction: the derection of the plateau end point 
    % ex) - 1: 100 W end point to start point during 100-200W trials -
    % Backmethod
    %     - 2: 200 W start point to end point during 100-200W tiral:
    %     Forward Method
    
    N = length(time);
    spm = round(N/abs(time(end)-time(1))); % the number of sampling per minute
    plateau_true = 0;
    s = time_window_min * spm; % the number of sampling per time window
        
    [~,m] = min(abs(time_start - time)); % finding a time stamp closes with the time start
    [~,n] = min(abs(time_end - time));   
    
    if direction == 2
        for i= m:(n-s)     
            a = polyfit(time(i:n), data_value(i:n), 1);
            if a(1) <= threshold && a(1) >= -1*threshold && peak2peak(data_value(m:i))<= 5 %0.1*abs(time(m)-time(i))
               plateau_true = 1; 
               plateau_start_index = i;
               plateau_end_index = n;
               break; 
            end
            if plateau_true == 0 && i == n-s
                plateau_start_index = i;
                plateau_end_index = n;
            end
        end
    end
    
    if direction == 1
        %direction
        for i = n:-1:(m+s)
            a = polyfit(time(m:i), data_value(m:i),1);
            if a(1) <= threshold && a(1) >= -1*threshold && peak2peak(data_value(m:i))<= 5 %0.1*abs(time(m)-time(i))
                plateau_true = 1; 
                plateau_start_index = m;
                plateau_end_index = i;
                break;
            end
            if plateau_true == 0 && i == m+s
                plateau_start_index = m;
                plateau_end_index = i;
            end
        end
    end
            
    if abs(time_start - time_end) < time_window_min
        
        plateau_start_index = nan;
        plateau_end_index = nan;        
    end
end