function [matches matches_entries] = aaaa_match_aeronet_with_modis(satellite)

% matches have the following format:
% [AERONET file, AERONET entry, MODIS file, MODIS entry]
% 
% matches_entries have the following format:
% [AERONET entry (40 columns), MODIS entry (21 columns)]
%
% NOT any more, now it has 16 entries, described in cell array
% fields_modis_matched

THRESHOLD_SPACE = 5; % kilometers
THRESHOLD_TIME = 30; % minutes

expected_matches_size = 2000;
found_matches = 0;
matches = zeros(expected_matches_size, 4);
matches_entries = zeros(expected_matches_size, 16);
% matches_entries = zeros(expected_matches_size, 61);

aeronet_files = dir('AERONET_data_from_modis_raw_processed\*.mat');
counter = 0;
for aeronet_site_index = 1 : length(aeronet_files)
    if (aeronet_files(aeronet_site_index).name(1) == '_')
        continue;
    end;
    
    % load the current locations coordinates to check if in wanted region
    load(['AERONET_data_from_modis_raw_processed\', aeronet_files(aeronet_site_index).name], 'longitude', 'latitude');
    if (~isInNorthernAmerica(longitude, latitude))
        continue;
    end;
    
    % this location is in our region, load the whole AERONET and MODIS
    % locations to find which ones to compare
    load(['AERONET_data_from_modis_raw_processed\', aeronet_files(aeronet_site_index).name], 'aeronet_data_better');
    load(['Satellite_MODIS_', satellite, '_data_raw_processed\051\aaa_per_location_mat_file\_fields_and_locs.mat'], 'locations');
    
    for loc = 1 : size(locations, 1)
        % for each MODIS locations check if they are close in space
        if (~areCloseInSpace(longitude, latitude, locations(loc, 1), locations(loc, 2), THRESHOLD_SPACE))
            continue;
        end;
        
        % now that we know that AERONET and MODIS are close in nearby, load
        % MODIS
        load(['Satellite_MODIS_', satellite, '_data_raw_processed\051\aaa_per_location_mat_file\file_num', num2str(loc), '.mat'], ['MODIS_', satellite, '_mat']);
        eval(['MODIS_data = MODIS_', satellite, '_mat;']);
        
        % now compare each entry and check if they are nigh
        match_found_at = 1;
        for i = 1 : size(aeronet_data_better, 1)
            current_aeronet = aeronet_data_better(i, :);
            
            for j = match_found_at : size(MODIS_data, 1)
                current_modis = MODIS_data(j, :);
                
                if (areCloseInTime(current_aeronet(1 : 5), current_modis(1 : 5), THRESHOLD_TIME))
                    % we have a match, save it
                    found_matches = found_matches + 1;
                    if (found_matches > length(matches))
                        matches = [matches; zeros(expected_matches_size, 4)];
                        matches_entries = [matches_entries; zeros(expected_matches_size, 16)];
%                         matches_entries = [matches_entries; zeros(expected_matches_size, 61)];
                    end;
                    
                    matches(found_matches, :) = [aeronet_site_index, i, loc, j];
                    matches_entries(found_matches, :) = [current_modis([1:7, 13, 16:18]), current_aeronet([35, 31, 25, 11, 9])];
%                     matches_entries(found_matches, :) = [current_aeronet, current_modis];
                    match_found_at = j;
                    break;
                elseif (isSecondTimeInFuture(current_aeronet(1 : 3), current_modis(1 : 3)))
                    break;
                end;
            end;
        end;
    end;
    counter = counter + 1;
end;
% counter       % this is number of sites in North America
matches = matches(1 : found_matches, :);
matches_entries = matches_entries(1 : found_matches, :);

% sort matches by AERONET entries
temps = matches(:, 1 : 2) * [10000000; 1];
[~, sorted] = sort(temps);
matches = matches(sorted, :);
matches_entries = matches_entries(sorted, :);

% now find duplicate entries, possible for various stupid reasons
duplicates = false(1, found_matches);
for i = 1 : found_matches - 1
    if (isSame(matches(i, 1 : 2), matches(i + 1, 1 : 2)) || isSame(matches(i, 3 : 4), matches(i + 1, 3 : 4)))
        duplicates(i + 1) = 1;
    end;
end;
matches = matches(~duplicates, :);
matches_entries = matches_entries(~duplicates, :);

% and again, for some additional stupid reasons
% now find duplicate entries, possible for various stupid reasons
duplicates = false(1, size(matches, 1));
for i = 1 : size(matches, 1) - 1
    if (isSame(matches(i, 1 : 2), matches(i + 1, 1 : 2)) || isSame(matches(i, 3 : 4), matches(i + 1, 3 : 4)))
        duplicates(i + 1) = 1;
    end;
end;
matches = matches(~duplicates, :);
matches_entries = matches_entries(~duplicates, :);

% and again, for some additional stupid reasons
% now find duplicate entries, possible for various stupid reasons
duplicates = false(1, size(matches, 1));
for i = 1 : size(matches, 1) - 1
    if (isSame(matches(i, 1 : 2), matches(i + 1, 1 : 2)) || isSame(matches(i, 3 : 4), matches(i + 1, 3 : 4)))
        duplicates(i + 1) = 1;
    end;
end;
matches = matches(~duplicates, :);
matches_entries = matches_entries(~duplicates, :);

function inNA = isInNorthernAmerica(long, lat)
% returns true (1) if the location is in North America
inNA = (((long < -50) & (long > -135)) & ((lat < 60) & (lat > 15)));

function areClose = areCloseInSpace(long1, lat1, long2, lat2, threshold)
% return (1) if locations are closer than THRESHOLD kilometers
areClose = (pos2dist(lat1, long1, lat2, long2, 1) <= threshold);

function areClose = areCloseInTime(time1, time2, threshold)
% return (1) if times are closer than THRESHOLD minutes, timeX variable is
% formatted in the following way:
% [year, month, day, hour, minute]
areClose = false;
if (sum(time1(1:3) == time2(1:3)) == 3)
    if (abs(time1(4 : 5) * [60; 1] - time2(4 : 5) * [60; 1]) < threshold)
        areClose = true;
    end;
end;

function inFuture = isSecondTimeInFuture(time1, time2)
% returns 1 if time2 is in future
inFuture = ((time1 * [366; 31; 1]) < (time2 * [366; 31; 1]));

function the_same = isSame(x, y)
the_same = (length(x) == sum((x == y) | (isnan(x) .* isnan(y))));