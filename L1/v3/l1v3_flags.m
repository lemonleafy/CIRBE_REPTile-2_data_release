%% Manually set the data flags for the v3 release files. 
% A complete list of dates that were flagged or corrected is available in
% invalid times.txt
clear
out_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main\CIRBE_REPTile-2_L1_v3';
flag1 = zeros(14704,1);
flag1(2567:14704)=1;

flag2 = ones(14076,1);
flag3 = ones(14260,1);
flag4 = ones(6896,1);
flag5 = ones(6628,1);
flag6 = ones(5107,1);

time7 = length(double(ncread(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230519v3_0.nc'),'Epoch')));
flag7 = zeros(time7,1);
flag7(1:1031)=1;

flag8 = zeros(6979,1);
flag8(463:6979) = 1;

flag9 = ones(2281,1);
flag10 = ones(29654,1);

flag11 = zeros(5214,1);
flag11(4491:5214) = 1;

ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230513v3_0.nc'),'invalid_data_flag',flag1)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230514v3_0.nc'),'invalid_data_flag',flag2)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230515v3_0.nc'),'invalid_data_flag',flag3)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230516v3_0.nc'),'invalid_data_flag',flag4)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230517v3_0.nc'),'invalid_data_flag',flag5)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230518v3_0.nc'),'invalid_data_flag',flag6)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230519v3_0.nc'),'invalid_data_flag',flag7)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230529v3_0.nc'),'invalid_data_flag',flag8)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230530v3_0.nc'),'invalid_data_flag',flag9)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230608v3_0.nc'),'invalid_data_flag',flag10)
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20230623v3_0.nc'),'invalid_data_flag',flag11)


flag12 = zeros(17165,1);
flag12(17165)=1;
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20240311v3_0.nc'),'invalid_data_flag',flag12)


flag13 = zeros(13163,1);
flag13(12543:13163)=1;
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20240217v3_0.nc'),'invalid_data_flag',flag13)


flag14 = zeros(26717,1);
flag14(14420:15527)=1;
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20240313v3_0.nc'),'invalid_data_flag',flag14)


flag15 = zeros(42355,1);
flag15(13480:13646)=1;
ncwrite(fullfile(out_dir,'CIRBE_REPTile-2_L1_20240802v3_0.nc'),'invalid_data_flag',flag15)