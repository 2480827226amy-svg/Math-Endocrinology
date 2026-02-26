%% Dosing function iteration 1
% Piecewise function
% Not to be incorporated into the model, just good for visualization

function y = tdosefunxx(t)
    Tdose = 3*10^6;
    y = zeros(size(t));  % Initialize output array
    
    for i = 0:52
        condition1 = (7*i < t) & (3+7*i > t);
        condition2 = (3+7*i <= t) & (5+7*i > t);
        condition3 = (5+7*i <= t) & (7+7*i >= t);
        
        % Apply piecewise conditions
        y(condition1) = Tdose;
        y(condition2) = 0.3 * Tdose;
        y(condition3) = 0;
    end
end

% Example usage and plot
t = linspace(0, 30);
tdose = tdosefunxx(t);

figure;
plot(t, tdose);
xlabel('t');
ylabel('Tdose');
title('Plot of tdosefunxx');
grid on;
