rng('shuffle')
%PATH = '//storage.neurologie.uni-tuebingen.de/bbnp_lab/Experimental Data/2023-01 HighDef/Conventional/Planning.xlsx';
PATH = 'C:/Users/BNP Lab/Documents/HighDef-temp/Planning.csv';
if isfile(PATH)
    T = readtable(PATH);
else
    T = [];
end


ID = input('Specify subject ID\n > ', 's');
%if ismember(ID, T.Subject)
%    for answer = input('This subject already has been planned! Overwrite? [y/n] ', 's')
%        if strcmpi(answer, 'y')
%            break
%        end
%    end
%end

% Select random order for left/right hotspot in initiation session:
InitiationConditions = cell(1,2);
r1 = round(rand); % Bernoulli 0 or 1 (very slight numeric bias towards 1; approx 0.4999 instead of 0.5)
InitiationConditions{1+r1} = 'L';
InitiationConditions{2-r1} = 'R';

InhibitionConditions = cell(1,2);
r2 = round(rand);
InhibitionConditions{1+r2} = 'L';
InhibitionConditions{2-r2} = 'R';

subjectRow = table({ID}, InitiationConditions(1), InitiationConditions(2), InhibitionConditions(1), InhibitionConditions(2), 'VariableNames', {'Subject', 'x1st_initiation', 'x2nd_initiation', 'x1st_bilateral', 'x2nd_bilateral'});
T = [T; subjectRow];

writetable(T, PATH);