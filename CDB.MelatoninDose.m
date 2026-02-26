%% Dr. Cecilia Diniz Behn's work on introducing exogenous melatonin
% Possibly helpful for a T dosing function

function dose_ex_mel = CDBmodel(t, melatonin_timing, melatonin_dosage)
    if isempty(melatonin_timing)
        dose_ex_mel = 0; % set exogenous melatonin to zero
        return;
end
    
    if t > 0 && t < 24
        t = mod(t, 24);
        if melatonin_timing - 0.1 <= t && t <= melatonin_timing + 0.3
            sigma = sqrt(0.002);
            mu = melatonin_timing + 0.1;
            
            ex_mel = (1 / (sigma * sqrt(2 * pi))) * exp(-(t - mu)^2 / (2 * sigma^2)); % Gaussian function
            
            x = 0:0.01:1;
            melatonin_values = max(x, sigma); % Assuming max_value is defined similarly in MATLAB
            max_value = max(melatonin_values);
            
            converted_dose = mg_conversion(melatonin_dosage); % Assuming mg_conversion is defined similarly in MATLAB
            
            normalize_ex_mel = (1 / max_value) * ex_mel; % normalize the values so the max is 1
            dose_ex_mel = converted_dose * normalize_ex_mel; % multiply by the dosage so the max = dosage
            
            return;
        else
            dose_ex_mel = 0;
            return;
        end
    else
        dose_ex_mel = 0;
        return;
    end
end

function melatonin_values = max(x, sigma)
    % Placeholder function for max_value
    melatonin_values = (1 / (sigma * sqrt(2 * pi))) * exp(-(x - 0.5).^2 / (2 * sigma^2)); % Example Gaussian values
end

function converted_dose = mg_conversion(melatonin_dosage)
    % Placeholder function for mg_conversion
    converted_dose = melatonin_dosage * 1; % Example conversion
end
