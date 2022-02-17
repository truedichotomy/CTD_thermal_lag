

%%
% pair_id = 130; % bad
% pair_id = 274;% good
% pair_id = 585;% bad
% pair_id = 1089; % good
% pair_id = 2120; % bad

plot(downcast(pair_id).date_num, -downcast(pair_id).depth, 'r-')
hold on
plot(downcast(pair_id+1).date_num, -downcast(pair_id+1).depth,'r--')
hold on
plot(upcast(pair_id).date_num, -upcast(pair_id).depth,'b-')
hold on
plot(upcast(pair_id+1).date_num, -upcast(pair_id+1).depth,'b--')
datetick
ylabel('z (m)')
xlabel('time')
title({'example of problematic profile identification', ['pair id = ' num2str(pair_id)]})
set(gca,'linewidth',1, 'fontsize', 14)

%%
for pair_id = 1089:1113 % good

plot(downcast(pair_id).date_num, -downcast(pair_id).depth, 'r.')
hold on
plot(upcast(pair_id).date_num, -upcast(pair_id).depth,'b.')
hold on
end

datetick
ylabel('z (m)')
xlabel('time')
set(gca,'linewidth',1, 'fontsize', 14)
