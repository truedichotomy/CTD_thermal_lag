%% load data
% load sylvia-20180501T0000-profile-sci-delayed.mat
% 
% glider_sci = sylvia_20180501T0000_profile_sc;
% clear sylvia_20180501T0000_profile_sc

%%
%
% select good (not a nan) data

% % time stamps of conductivity are same as those of tempature
% good_id0 = find(~isnan(glider_sci.conductivity));
good_id1 = find(~isnan(glider_sci.temperature));
[sci_data.time, good_id2, id3] = unique(glider_sci.time,'sorted');
good_id = intersect(good_id1, good_id2);

sci_data = IndexedStructCopy(glider_sci, good_id);

%%
% % raw time is in "seconds since 1970-01-01T00:00:00Z";
% smart way to convert time to date number format
sci_data.date_num = datenum(1970, 1, 1, 0, 0, sci_data.time); 

sci_data.z = -sci_data.depth;


%% definie function to extract data into desired subsets

clear segment_data
% need to find a better method to select segment
% selected by eyeball inspection (beginning and end of a segment)
segment_t1 = 1525769645.76123; % sci_data.time(53023)
segment_t2 = 1525780477.42612; % sci_data.time(57992)

segment_data = IndexedStructCopy(sci_data, sci_data.time>=segment_t1 & sci_data.time<=segment_t2);

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