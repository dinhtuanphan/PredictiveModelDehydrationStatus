function [time_average_s, time_average_e, average_s_to_e, average_last_1min, average_plateau_first_10min] = A_average(t_start,t_end, spm, time, data, i_stable, is_it_C)

    if ~isnan(is_it_C)
       t_start = sort([t_start time(i_stable)]);
       m = sum(t_start <= time(i_stable));
       t_start = t_start(m:end);
       t_end = t_end(m-1:end);
    end
      
    time_average_s = t_start;
    time_average_e = t_end;
    
       
    N = length(t_start);    
    i_start = nan(N,1);
    i_end = nan(N,1); 
            
    for i = 1:N
                [~, i_start(i)] = min(abs(t_start(i) - time));
        [~, i_end(i)] = min(abs(t_end(i)-time));
    end
    
    average_s_to_e = nan(N, 1);
    average_last_1min = nan(N, 1);
    
    for i = 1:N
        average_s_to_e(i) = mean(data(i_start(i):i_end(i))); 
        average_plateau_first_10min = mean(data(i_stable:i_stable+spm*10));
        average_last_1min(i) = mean(data(i_end(i)-spm:i_end(i)));
    end

end
   