% Fix: range proton channels 1 and 2 were combined through 2023-07-30.
% Channel 1 recorded nothing (empty), channel 2 recorded the sum of both
% channels' counts. Set channel 1 to the standard fill value for every
% day up to and including 2023-07-30. Channel 2 is left untouched but
% is documented in the release notes as containing combined
% ch1+ch2 counts for this period.

clear

out_dir = 'C:\Users\wesle\Documents\LASP research\Scripts\CIRBE_REPTile-2_data_release-main\CIRBE_REPTile-2_L1_v3';

fill_value = single(9.9692100e+36);   % matches Pbins_RNG FillValue in nccreate

trange = ['2023-04-19'; '2023-07-30'];   % inclusive cutoff
start_datenum = datenum(trange(1,:));
end_datenum   = datenum(trange(2,:));

n_fixed = 0;
n_skipped = 0;

for i = start_datenum:end_datenum

    filedate = datestr(i, 'yyyymmdd');
    target_file = fullfile(out_dir, ['CIRBE_REPTile-2_L1_', filedate, 'v3_0.nc']);

    if exist(target_file, 'file') ~= 2
        n_skipped = n_skipped + 1;
        continue
    end

    Pbins_RNG = ncread(target_file, 'Pbins_RNG');   % [50, num_t]

    if size(Pbins_RNG,1) ~= 50
        error('%s: Pbins_RNG has %d channels, expected 50.', filedate, size(Pbins_RNG,1))
    end

    Pbins_RNG(1,:) = fill_value;

    ncwrite(target_file, 'Pbins_RNG', Pbins_RNG)

    % Verify
    check = ncread(target_file, 'Pbins_RNG');
    if ~all(isnan(check(1,:)))
        % ncread converts values matching the variable's _FillValue
        % attribute to NaN, so after write-back check should read as NaN.
        error('%s: channel 1 fill verification failed.', filedate)
    end

    fprintf('Fixed: %s (channel 1 set to fill value)\n', filedate)
    n_fixed = n_fixed + 1;
end

fprintf('\nDone. Fixed %d files, skipped %d (no file present).\n', n_fixed, n_skipped)