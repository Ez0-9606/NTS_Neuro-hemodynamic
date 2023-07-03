function [occu, info] = get_information(state_data, raw_fr)
    % mean firing rate for each state and each neuron
    fr = zeros(max(state_data) - min(state_data) + 1, size(raw_fr,2));
    % occupancy ratio
    occu = zeros(max(state_data) - min(state_data) + 1, 1);
    
    for i = 1:numel(occu)
        temp_idx = find(state_data == (min(state_data)+i-1));
        occu(i) = numel(temp_idx)/numel(state_data);
        if isempty(temp_idx)
            fr(i,:) = 0;
        else
%             fr(i,:) = mean(round(20*raw_fr(temp_idx,:)), 1);
            fr(i,:) = mean(round(raw_fr(temp_idx,:)), 1);
        end
    end
    
    % mean firing rate across states for each neuron
    avg_fr = sum(fr.*repmat(occu,1, size(raw_fr,2)), 1);
    
    % Information content for ecah neuron
    info = zeros(size(raw_fr,2), 1);
    for i = 1:numel(info)
        temp = 0;
        for j = 1:numel(occu)
            if fr(j,i)*occu(j) == 0
                continue;
            else
                val = (fr(j,i)/avg_fr(i))*log2(fr(j,i)/avg_fr(i))*occu(j);
            end
%             if isnan(val)|| isinf(val)
            temp = temp + val;
        end
        info(i) = temp;
    end
end