% Fix: indexing dimensionality bug in v3 proton AND electron data for the
% 4 special days. Original v2 production used linear indexing
% (Var(indexes)) on [50,num_t]/[10,num_t] arrays, which does not select
% whole time-columns. Correct form is Var(:, indexes). This corrects all
% four count arrays (Pbins_RNG, Pbins_PEN, Ebins_RNG, Ebins_PEN) by
% reloading raw counts from source and reindexing properly.

clear

out_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main\CIRBE_REPTile-2_L1_v3';
sp_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main\special indexes';
internal_dir = 'Z:\combined_science';

special_days = struct( ...
    'date', {'2023_06_21', '2023_06_25', '2023_07_01', '2023_11_05'}, ...
    'filedate', {'20230621', '20230625', '20230701', '20231105'}, ...
    'indexfile', { ...
        fullfile(sp_dir, 'special_indexes_2023_06_21.mat'), ...
        fullfile(sp_dir, 'special_indexes_2023_06_25.mat'), ...
        fullfile(sp_dir, 'special_indexes_2023_07_01.mat'), ...
        fullfile(sp_dir, 'special_indexes_2023_11_05.mat')});

for s = 1:numel(special_days)
    ds = special_days(s).date;
    fd = special_days(s).filedate;

    target_file = fullfile(out_dir, ['CIRBE_REPTile-2_L1_', fd, 'v3_0.nc']);
    source_file = fullfile(internal_dir, ['CIRBE_L1_combined_science_', ds, '_V0.nc']);

    if exist(target_file, 'file') ~= 2
        error('Target file not found: %s', target_file)
    end
    if exist(source_file, 'file') ~= 2
        error('Source file not found: %s', source_file)
    end

    load(special_days(s).indexfile, 'indexes');

    num_t = length(ncread(target_file, 'Epoch'));

    Pbins_RNG_raw = ncread(source_file, 'pcounts');
    Pbins_PEN_raw = ncread(source_file, 'hpcounts');
    Ebins_RNG_raw = ncread(source_file, 'Ecounts');
    Ebins_PEN_raw = ncread(source_file, 'hEcounts');

    if size(Pbins_RNG_raw,1) ~= 50 || size(Pbins_PEN_raw,1) ~= 10 || ...
       size(Ebins_RNG_raw,1) ~= 50 || size(Ebins_PEN_raw,1) ~= 10
        error(['%s: unexpected source array shape — pcounts %dx%d, hpcounts %dx%d, ', ...
               'Ecounts %dx%d, hEcounts %dx%d.'], ds, ...
            size(Pbins_RNG_raw,1), size(Pbins_RNG_raw,2), ...
            size(Pbins_PEN_raw,1), size(Pbins_PEN_raw,2), ...
            size(Ebins_RNG_raw,1), size(Ebins_RNG_raw,2), ...
            size(Ebins_PEN_raw,1), size(Ebins_PEN_raw,2))
    end

    Pbins_RNG_fixed = Pbins_RNG_raw(:, indexes);
    Pbins_PEN_fixed = Pbins_PEN_raw(:, indexes);
    Ebins_RNG_fixed = Ebins_RNG_raw(:, indexes);
    Ebins_PEN_fixed = Ebins_PEN_raw(:, indexes);

    fixed = {Pbins_RNG_fixed, Pbins_PEN_fixed, Ebins_RNG_fixed, Ebins_PEN_fixed};
    names = {'Pbins_RNG','Pbins_PEN','Ebins_RNG','Ebins_PEN'};
    for k = 1:numel(fixed)
        if size(fixed{k},2) ~= num_t
            error('%s: reindexed %s has %d time points; Epoch has %d.', ...
                ds, names{k}, size(fixed{k},2), num_t)
        end
    end

    ncwrite(target_file, 'Pbins_RNG', Pbins_RNG_fixed)
    ncwrite(target_file, 'Pbins_PEN', Pbins_PEN_fixed)
    ncwrite(target_file, 'Ebins_RNG', Ebins_RNG_fixed)
    ncwrite(target_file, 'Ebins_PEN', Ebins_PEN_fixed)

    % Verify writes
    check_Pbins_RNG = ncread(target_file, 'Pbins_RNG');
    check_Pbins_PEN = ncread(target_file, 'Pbins_PEN');
    check_Ebins_RNG = ncread(target_file, 'Ebins_RNG');
    check_Ebins_PEN = ncread(target_file, 'Ebins_PEN');

    if ~isequal(check_Pbins_RNG, Pbins_RNG_fixed) || ...
       ~isequal(check_Pbins_PEN, Pbins_PEN_fixed) || ...
       ~isequal(check_Ebins_RNG, Ebins_RNG_fixed) || ...
       ~isequal(check_Ebins_PEN, Ebins_PEN_fixed)
        error('%s: write verification failed after reindexing fix.', ds)
    end

    fprintf('Fixed and verified: %s (protons + electrons)\n', ds)

    clear indexes
end

disp('All 4 special-day v3 files reindexed correctly for both species.')