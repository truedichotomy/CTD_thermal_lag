%% load data
% load ru28-20170424T1310-profile-sci-delayed_0823_17a1_4ff7.mat
% glider_sci = ru28_20170424T1310_profile_sci_;
% clear ru28_20170424T1310_profile_sci_

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


%%

% assuming 0.5 Hz sampling frequency, 
% this moving mean window corresponds to 10 seconds time window
% how is this different from low-pass filtering? MERCKELBACH 2021 paper

sci_data.pressure_smooth = smoothdata(sci_data.pressure, 'movmean', 5); 
sci_data.z_smooth = smoothdata(sci_data.z, 'movmean', 5); 


%%
% profile id does not increase monotonically. try to use profile time as an index to differentiate casts.
% plot(glider_sci.time, glider_sci.profile_id)
% plot(glider_sci.time, glider_sci.profile_time)


%% definie function to extract data into desired subsets
segment_t1 = 1494178330.85495;
segment_t2 = 1494185750.71231;

segment_data = IndexedStructCopy(sci_data, sci_data.time>=segment_t1 & sci_data.time<=segment_t2);

%%
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