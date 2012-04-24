function [matches_aeronet_terra_misr, found_misr_matches] = aaaa_match_aeronet_with_both_terra_and_misr(matches_terra, matches_misr)

% matches have the following format:
% [AERONET file, AERONET entry, MODIS file, MODIS entry, MISR file, MISR entry]

matches_aeronet_terra_misr = zeros(min([size(matches_terra, 1), size(matches_misr, 1)]), 6);
match_counter = 0;
found_misr_matches = zeros(1, size(matches_misr, 1));

last_found_at = 1;
for i = 1 : size(matches_misr, 1)
    current_misr_match = matches_misr(i, :);
    for j = last_found_at : size(matches_terra, 1)
        current_terra_match = matches_terra(j, :);
        
        % if they have the same match in AERONET entries
        if (isSame(current_misr_match(1 : 2), current_terra_match(1 : 2)))
            match_counter = match_counter + 1;
            matches_aeronet_terra_misr(match_counter, :) = [current_terra_match, current_misr_match(3 : 4)];
            
            % move last_found_at only if MISR entry found for the first 
            % time; WHY you ask? since there might be more than 1 matches
            if (found_misr_matches(i) == 0)
                last_found_at = j;
            end;
            
            found_misr_matches(i) = 1;
        else
            if (found_misr_matches(i))
                break;
            end;
        end;
    end;
    
    if (mod(i, 100) == 0)
        i
    end;
end;
matches_aeronet_terra_misr = matches_aeronet_terra_misr(1 : match_counter, :);

function the_same = isSame(x, y)
the_same = (length(x) == sum((x == y)));
% the_same = (length(x) == sum((x == y) | (isnan(x) .* isnan(y))));