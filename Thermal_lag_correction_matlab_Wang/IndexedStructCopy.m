%% definie function to extract data into desired subsets
% function from https://www.mathworks.com/matlabcentral/answers/405944-how-do-i-extract-subset-of-all-fields-from-structure
function T = IndexedStructCopy(S, Condition, FieldList)
    if nargin == 2
        FieldList = fieldnames(S);
    end
    for iField = 1:numel(FieldList)
        Field    = FieldList{iField};
        T.(Field) = S.(Field)(Condition);
    end
end