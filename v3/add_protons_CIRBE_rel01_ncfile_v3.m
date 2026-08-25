%% Copy the v2 release files, add proton counts from the lasp store files, save new files as v3
clear
tic
src_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main';
out_dir = fullfile(src_dir, 'CIRBE_REPTile-2_L1_v3');

trange = ['2023-04-19';'2024-09-28'];  % full public v2 range
start_datenum = datenum(trange(1,:));
end_datenum = datenum(trange(2,:));

sp_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main\special indexes\';
special_days = struct( ... 
    'date', {'2023_06_21', '2023_06_25', '2023_07_01', '2023_11_05'}, ...
    'indexfile', {[sp_dir, 'special_indexes_2023_06_21.mat'], [sp_dir, 'special_indexes_2023_06_25.mat'], ...
    [sp_dir, 'special_indexes_2023_07_01.mat'], [sp_dir, 'special_indexes_2023_11_05.mat']} ...
); % Indexes selecting the usable data from four partially corrupted days

for i = start_datenum:end_datenum
    file_datestring = datestr(i,"yyyy_mm_dd");
    filedate = datestr(i,"yyyymmdd");

    file_path_1 = ['Z:\combined_science\CIRBE_L1_combined_science_',file_datestring,'_V0.nc'];  % internal source
    old_file = fullfile(src_dir, ['CIRBE_REPTile-2_L1_v2/CIRBE_REPTile-2_L1_',filedate,'v2_0.nc']);
    new_file = fullfile(out_dir, ['CIRBE_REPTile-2_L1_',filedate,'v3_0.nc']);

    if exist(file_path_1,'file') ~= 2 || exist(old_file,'file') ~= 2
        continue
    end

    disp(file_datestring)

    copyfile(old_file, new_file)   % never touch old_file itself

    Pbins_RNG = ncread(file_path_1,"pcounts");
    Pbins_PEN = ncread(file_path_1,"hpcounts");

    is_special = false;
    for s = 1:numel(special_days)
        if strcmp(special_days(s).date, file_datestring)
            is_special = true;
            load(special_days(s).indexfile, 'indexes');
            Pbins_RNG = Pbins_RNG(:, indexes);
            Pbins_PEN = Pbins_PEN(:, indexes);
            break
        end
    end

    published_num_t = length(ncread(new_file,"Epoch"));
    if size(Pbins_RNG,2) ~= published_num_t
        error(['Length mismatch on ', file_datestring, ...
               ': proton source has ', num2str(size(Pbins_RNG,2)), ...
               ' points (special=', num2str(is_special), '), published Epoch has ', num2str(published_num_t)])
    end

    ncwrite(new_file,'Pbins_RNG',Pbins_RNG)
    ncwrite(new_file,'Pbins_PEN',Pbins_PEN)

    ncwriteatt(new_file,"/","Data_version",'3.0')
end
toc