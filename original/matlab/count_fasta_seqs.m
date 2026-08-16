cd E:\FASTA\Gabriel\

% Sample data
data = readtable('count_seqs_fasta.csv', 'Delimiter', ',');

% Add new columns to the existing table
data.StartDay = zeros(size(data, 1), 1);
data.StartMonth = zeros(size(data, 1), 1);
data.StartYear = zeros(size(data, 1), 1); 
data.EndDay = zeros(size(data, 1), 1);
data.EndMonth = zeros(size(data, 1), 1);
data.EndYear = zeros(size(data, 1), 1);

%%
% Loop through each row of the data
for i = 1:size(data, 1)
    % Extract start and end date strings
    dateStr = char(data.FileName(i));
    startMatch = regexp(dateStr, 'Start_(\d+-[a-zA-Z]+-\d+)_', 'tokens', 'once');
    endMatch = regexp(dateStr, 'End_(\d+-[a-zA-Z]+-\d+)_', 'tokens', 'once');
    
     % Check if matches are found
    if ~isempty(startMatch) && ~isempty(endMatch)
        % Extract start and end date strings from the matches
        startStr = startMatch{1};
        endStr = endMatch{1};
        
        % Convert start and end date strings to datetime objects
        startDate = datetime(startStr, 'InputFormat', 'dd-MMM-yyyy');
        endDate = datetime(endStr, 'InputFormat', 'dd-MMM-yyyy');
        
        % Update the table with the extracted information
        data.StartDay(i) = day(startDate);
        data.StartMonth(i) = month(startDate);
        data.StartYear(i) = year(startDate); 
        data.EndDay(i) = day(endDate);
        data.EndMonth(i) = month(endDate);
        data.EndYear(i) = year(endDate);
    else
        % Handle the case where no match is found
        disp(['Error: No match found for row ', num2str(i)]);
    end
end

