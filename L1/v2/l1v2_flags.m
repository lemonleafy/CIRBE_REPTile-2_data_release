flag1 = zeros(14704,1);
flag1(2567:14704)=1;

flag2 = ones(14076,1);
flag3 = ones(14260,1);
flag4 = ones(6896,1);
flag5 = ones(6628,1);
flag6 = ones(5107,1);

time7 = length(double(ncread("CIRBE_REPTile-2_L1_20230519v2_0.nc","Epoch")));
flag7 = zeros(time7,1);
flag7(1:1031)=1;

flag8 = zeros(6979,1);
flag8(463:6979) = 1;

flag9 = ones(2281,1);
flag10 = ones(29654,1);

flag11 = zeros(5214,1);
flag11(4491:5214) = 1;

ncwrite('CIRBE_REPTile-2_L1_20230513v2_0.nc','invalid_data_flag',flag1)
ncwrite('CIRBE_REPTile-2_L1_20230514v2_0.nc','invalid_data_flag',flag2)
ncwrite('CIRBE_REPTile-2_L1_20230515v2_0.nc','invalid_data_flag',flag3)
ncwrite('CIRBE_REPTile-2_L1_20230516v2_0.nc','invalid_data_flag',flag4)
ncwrite('CIRBE_REPTile-2_L1_20230517v2_0.nc','invalid_data_flag',flag5)
ncwrite('CIRBE_REPTile-2_L1_20230518v2_0.nc','invalid_data_flag',flag6)
ncwrite('CIRBE_REPTile-2_L1_20230519v2_0.nc','invalid_data_flag',flag7)
ncwrite('CIRBE_REPTile-2_L1_20230529v2_0.nc','invalid_data_flag',flag8)
ncwrite('CIRBE_REPTile-2_L1_20230530v2_0.nc','invalid_data_flag',flag9)
ncwrite('CIRBE_REPTile-2_L1_20230608v2_0.nc','invalid_data_flag',flag10)
ncwrite('CIRBE_REPTile-2_L1_20230623v2_0.nc','invalid_data_flag',flag11)