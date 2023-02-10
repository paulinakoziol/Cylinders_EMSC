function [pos_1, pos_2] = find_value_min_max(vector, min_value, max_value)

[output_value1, pos_1] = min(abs(vector-min_value));

[output_value2, pos_2] = min(abs(vector-max_value));


end % End of find_value_min_max



