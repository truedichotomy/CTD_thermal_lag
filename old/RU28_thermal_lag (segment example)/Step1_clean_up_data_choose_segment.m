%% load data
% load ru28-20170424T1310-profile-sci-delayed_0823_17a1_4ff7.mat
% glider_sci = ru28_20170424T1310_profile_sci_;
% clear ru28_20170424T1310_profile_sci_


load RU28_glider_science_data_raw.mat

% convert source file names from char to string, in preparation for data extraction
glider_sci.source_file = string(glider_sci.source_file);
glider_sci.trajectory = string(glider_sci.trajectory);

%%
%
% select good (not a nan) data

% % time stamps of conductivity are same as those of tempature
% good_id0 = find(~isnan(glider_sci.conductivity));
% good_id0 = find(~isnan(glider_sci.temperature));

good_id1 = find(~isnan(glider_sci.sci_ctd41cp_timestamp));
[~, good_id2, ~] = unique(glider_sci.sci_ctd41cp_timestamp,'sorted');
good_id = intersect(good_id1, good_id2);

% [glider_sci.time, good_id2, id3] = unique(glider_sci.time,'sorted');

sci_data = IndexedStructCopy(glider_sci, good_id);

sci_data = IndexedStructCopy(glider_sci, good_id);

sci_data.ctd_time = sci_data.sci_ctd41cp_timestamp;

% % raw time is in "seconds since 1970-01-01T00:00:00Z";
% smart way to convert time to date number format
sci_data.date_num = datenum(1970, 1, 1, 0, 0, sci_data.ctd_time); 

sci_data.z = -sci_data.depth;

%% extract data for each segment

% each unique source file name is asociated with a survey segment (between two surfacing points)
source_file_name = unique(sci_data.source_file,'rows');

clear segment

n_segment = size(source_file_name,1);
for iter = 1:n_segment
    segment(iter) = IndexedStructCopy(sci_data, ...
        sci_data.source_file == source_file_name(iter));
end

%%
segment_file = "ru28-2017-113-3-302-dbd(02220302)";
segment_data = IndexedStructCopy(sci_data, sci_data.source_file==segment_file);

% % old selection for demonstration
% segment_t1 = 1494178330.85495;
% segment_t2 = 1494185750.71231;
% segment_data = IndexedStructCopy(sci_data, sci_data.time>=segment_t1 & sci_data.time<=segment_t2);


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