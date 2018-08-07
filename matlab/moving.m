function y = moving(dat, n)
% Allocate space for vector y containing median values.
y = NaN(size(dat));

itr = length(dat);
for i = 1:itr
    % Determine bounds of window
    lb = i-n; ub = i+n;
    
    % Shrink window near edges of dataset
    if lb <= 0
        lb = 1;

        % Make size of window odd if it is even
        if mod(length(lb:ub),2) == 0
           ub = ub+1;
        end
    
    elseif ub > itr
        ub = itr;
        
        % Make size of window odd if it is even
        if mod(length(lb:ub),2) == 0
           ub = ub-1;
        end
    end
    
    % Calculate median of window
    y(i) = median(dat(lb:ub));
    clear lb ub
end