load("mode_times.mat") 
%trange = ['2024-01-01';'2024-04-15'];
trange = ['2023-06-08';'2023-06-08'];
%trange = ['2024-03-01';'2024-03-31'];
start_datenum = datenum(trange(1,:));
end_datenum = datenum(trange(2,:));

for i = start_datenum:end_datenum
    file_datestring = datestr(i,"yyyy_mm_dd");
    file_path_1 = ['CIRBE_L1_combined_science_',file_datestring,'_V0.nc'];
    file_path_2 = ['CIRBE_L1_Nom_science_',file_datestring,'_V0.nc'];
    fileID_CIRBE_1 = fopen(file_path_1);
    if fileID_CIRBE_1 > 0
        disp(file_datestring)
        filedate = datestr(i,"yyyymmdd");
        Datafilename = strjoin(["C:\Users\2001d\OneDrive\Desktop\LASP Data\Data release\netcdf file\daily_netcdf\CIRBE_REPTile-2_L1_",filedate,"v1_1.nc"],"");

        Alt0 = ncread(file_path_1,"Altitude");
        Lat0 = ncread(file_path_1,"Latitude");
        Lon0 = ncread(file_path_1,"Longitude");
        timeraw = double(ncread(file_path_2,"time_ymdhms"))';
        time = datetime(timeraw);
        %t = double(convertTo(time,'epochtime','Epoch','2000-01-01')); %seconds since 2000-01-01
        t = double(convertTo(time,'epochtime','Epoch','2000-01-01','TicksPerSecond',1000))/1000; %seconds since 2000-01-01

        kext = 0; % for L and MLT
        options = [0,0,0,0,0];
        sysaxes = 0;
        maginput = zeros(1,25);
        [L0_1,~,B,Bmin,~,MLT0] = onera_desp_lib_make_lstar(kext,options,sysaxes,time,Alt0,Lat0,Lon0,maginput);
        L0=single(abs(L0_1));

        HK_time_raw = double(ncread(file_path_1,"Pointing_timestamp")); % no ymdhms for offangle
        HK_time = datenum('1999/12/31 23:59:23')+(HK_time_raw/1000/3600/24);
        HK_time = datetime(HK_time,'ConvertFrom','datenum');
        %thk = double(convertTo(HK_time,'epochtime','Epoch','2000-01-01')); %seconds since 2000-01-01
        thk = double(convertTo(HK_time,'epochtime','Epoch','2000-01-01','TicksPerSecond',1000))/1000; %seconds since 2000-01-01

        num_t = length(time);
        num_tHK = length(HK_time);

        %t = zeros(num_t,1);
        %thk = zeros(num_tHK,1);

        Ebins_RNG = ncread(file_path_1,"Ecounts");
        Pbins_RNG = ncread(file_path_1,"pcounts");
        Ebins_PEN = ncread(file_path_1,"hEcounts");
        Pbins_PEN = ncread(file_path_1,"hpcounts");

        OffAng0 = single(ncread(file_path_1,"offpointing_mag_model"));

        D1 = ncread(file_path_1,"D1");
        D2 = ncread(file_path_1,"D2");
        D3 = ncread(file_path_1,"D3");
        D4 = ncread(file_path_1,"D4");
        D12n = ncread(file_path_1,"D12n");
        D123n = ncread(file_path_1,"D123n");
        D1234n = ncread(file_path_1,"D1234n");
        D1234 = ncread(file_path_1,"D1234");
        G = ncread(file_path_1,"G");
        IntePrd = ncread(file_path_1,"sci_intg_prd");

        xGEO = onera_desp_lib_rotate([Alt0(:), Lat0(:), Lon0(:)],'gdz2geo');
        xSM = onera_desp_lib_rotate(xGEO,'geo2sm',time);

        [azimuth,elevation,r] = cart2sph(xSM(:,1),xSM(:,2),xSM(:,3));

        mlt=(azimuth+pi)*12/pi;

        mlt(mlt<0)=mlt(mlt<0)+24;

        mlat=elevation*180/pi;

        valflag = zeros(num_t,1);

        spsn = double(convertTo(mode_times,'epochtime','Epoch','2000-01-01'));
        spm = t;
        for j = 1:num_t
            if ((spm(j) > spsn(1) & spm(j) < spsn(2)) | (spm(j) > spsn(3) & spm(j) < spsn(4)) | (spm(j) > spsn(5) & spm(j) < spsn(6)) | (spm(j) > spsn(7) & spm(j) < spsn(8)) | (spm(j) > spsn(9) & spm(j) < spsn(10)))
                spm(j) = 1;
            elseif(spm(j) > spsn(11) & spm(j) < spsn(12))
                spm(j) = 2;
            else
                spm(j) = 0;
            end
        end

        point_mode = single(spm);

        FileDir=dir(Datafilename);
        if isempty(FileDir)
        else
            delete ([Datafilename])
        end
        nccreate(Datafilename,'Epoch','Dimensions',{'Epoch',num_t},...
            'FillValue',-9223372036854775806,"Datatype","double")
        nccreate(Datafilename,'Alt','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'Lat','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'Lon','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'L','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'MLT','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'MLAT','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        % nccreate(Datafilename,'ChNum_RanP','Dimensions',{'ChNum_RanP',50},...
        %     'FillValue','disable')
        nccreate(Datafilename,'Pbins_RNG','Dimensions',{'ChNum_RanP',50,'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'Pbins_PEN','Dimensions',{'ChNum_PenP',10,'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        % nccreate(Datafilename,'ChNum_RanE','Dimensions',{'ChNum_RanE',50},...
        %     'FillValue',disable')
        nccreate(Datafilename,'Ebins_RNG','Dimensions',{'ChNum_RanE',50,'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'Ebins_PEN','Dimensions',{'ChNum_PenE',10,'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")

        nccreate(Datafilename,'D1','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D2','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D3','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D4','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'G','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D12n','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D123n','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D1234n','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'D1234','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'IntePrd','Dimensions',{'Epoch',num_t},...
            'FillValue',9.9692100e+36,"Datatype","single")
        nccreate(Datafilename,'OffAng_timestamp','Dimensions',{'OffAng_Ts',num_tHK},...
            'FillValue',-9223372036854775806,"Datatype","double")
        nccreate(Datafilename,'OffAng','Dimensions',{'OffAng_Ts',num_tHK},...
            'FillValue',9.9692100e+36,"Datatype","single")
        
        nccreate(Datafilename,'spacecraft_pointing_mode','Dimensions',{'Epoch',num_t},...
            'FillValue',999,"Datatype","single")
        nccreate(Datafilename,'invalid_data_flag','Dimensions',{'Epoch',num_t},...
            'FillValue',999,"Datatype","single")

        ncwrite(Datafilename,'Epoch',t)
        ncwrite(Datafilename,'Alt',Alt0)
        ncwrite(Datafilename,'Lat',Lat0)
        ncwrite(Datafilename,'Lon',Lon0)
        ncwrite(Datafilename,'L',L0)
        ncwrite(Datafilename,'MLT',mlt)
        ncwrite(Datafilename,'MLAT',mlat)

        ncwrite(Datafilename,'Ebins_RNG',Ebins_RNG)
        ncwrite(Datafilename,'Pbins_RNG',Pbins_RNG)
        ncwrite(Datafilename,'Ebins_PEN',Ebins_PEN)
        ncwrite(Datafilename,'Pbins_PEN',Pbins_PEN)
        ncwrite(Datafilename,'D1',D1)
        ncwrite(Datafilename,'D2',D2)
        ncwrite(Datafilename,'D3',D3)
        ncwrite(Datafilename,'D4',D4)
        ncwrite(Datafilename,'G',G)
        ncwrite(Datafilename,'D12n',D12n)
        ncwrite(Datafilename,'D123n',D123n)
        ncwrite(Datafilename,'D1234n',D1234n)
        ncwrite(Datafilename,'D1234',D1234)
        ncwrite(Datafilename,'IntePrd',IntePrd)

        ncwrite(Datafilename,'OffAng_timestamp',thk)
        ncwrite(Datafilename,'OffAng',OffAng0)
        ncwrite(Datafilename,"spacecraft_pointing_mode",point_mode)
        ncwrite(Datafilename,"invalid_data_flag",valflag)

    end
end


%%

%trange = ['2024-01-01';'2024-04-15'];
trange = ['2023-06-08';'2023-06-08'];
%trange = ['2024-03-01';'2024-03-31'];
start_datenum = datenum(trange(1,:));
end_datenum = datenum(trange(2,:));

varnames = string({'Epoch', 'Alt', 'Lat', 'Lon', 'L', 'MLT', 'MLAT', 'Pbins_RNG', 'Pbins_PEN', 'Ebins_RNG', 'Ebins_PEN', 'D1', 'D2', 'D3', 'D4', 'G', 'D12n', 'D123n', 'D1234n', 'D1234', 'IntePrd', 'OffAng_timestamp', 'OffAng', 'spacecraft_pointing_mode', 'invalid_data_flag'});

for i = start_datenum:end_datenum
    file_datestring = datestr(i,"yyyy_mm_dd");
    file_path_1 = ['CIRBE_L1_combined_science_',file_datestring,'_V0.nc'];
    file_path_2 = ['CIRBE_L1_Nom_science_',file_datestring,'_V0.nc'];
    fileID_CIRBE_1 = fopen(file_path_1);
    if fileID_CIRBE_1 > 0
        disp(file_datestring)
        filedate = datestr(i,"yyyymmdd");
        Datafilename = strjoin(["C:\Users\2001d\OneDrive\Desktop\LASP Data\Data release\netcdf file\daily_netcdf\CIRBE_REPTile-2_L1_",filedate,"v1_1.nc"],"");
        timeraw = double(ncread(file_path_2,"time_ymdhms"))';
        time = datetime(timeraw);
        t = int64(convertTo(time,'epochtime','Epoch','2000-01-01'));
        HK_time_raw = double(ncread(file_path_1,"Pointing_timestamp")); % no ymdhms for offangle
        HK_time = datenum('1999/12/31 23:59:23')+(HK_time_raw/1000/3600/24);
        HK_time = datetime(HK_time,'ConvertFrom','datenum');
        thk = int64(convertTo(HK_time,'epochtime','Epoch','2000-01-01'));
        num_t = length(time);
        num_tHK = length(HK_time);
        
        % universal attributes
        ncwriteatt(Datafilename,"/","Project",'NSF 3U Cubesat');
        ncwriteatt(Datafilename,"/",'Source_name','CIRBE>Colorado Inner Radiation Belt Experiment');
        ncwriteatt(Datafilename,"/",'Disicpline','Solar Physics>Heliospheric Physics, Space Physics>Magnetospheric Science');
        ncwriteatt(Datafilename,"/",'Descriptor','REPTile-2>Relativistic Electron and Proton Telescope Integrated Little Experiment-2');
        %ncwriteatt(Datafilename,"/",'File_naming_convention','source_datatype_descriptor');
        ncwriteatt(Datafilename,"/",'Data_version','1.0');
        ncwriteatt(Datafilename,"/",'Instrument_type','Particles (space)');
        ncwriteatt(Datafilename,"/",'Mission_group','Cubesats');
        ncwriteatt(Datafilename,"/",'PI_name','Xinlin Li');
        ncwriteatt(Datafilename,"/",'PI_affiliation','University of Colorado at Boulder');
        ncwriteatt(Datafilename,"/",'Acknowledgements','The REPTile-2 data are provided by the University of Colorado; the REPTile-2 PI  is Dr. Xinlin Li. The references in metadata REF1 and REF2 can be cited for the REPTile-2 data products. We acknowledge the use of the IRBEM library (4.4.0), the latest version of which can be found at https://doi.org/10.5281/zenodo.6867552.');
        ncwriteatt(Datafilename,"/",'Rules_of_use','The REPTile-2 data are available for science research. Users should contact the PI or other team member to discuss the appropriate use of these data. Users publishing results should provide appropriate acknowledgement, should provide the version number of the data being used, and could also offer co-authorship to PI or other team member.');
        ncwriteatt(Datafilename,"/",'TEXT','CIRBE is a 3U-CubeSat designed and developed by students and engineers at the Laboratory for Atmospheric and Space Physics. The primary objective of the science mission is to understand the formation of the inner belt (L<2) electrons (100s of keV to multiple MeV), and to determine the source, the intensity and dynamic variations of these electrons. The goal is to make accurate measurements with fine energy resolution (>40 channels) for electrons of 0.3-3 MeV throughout the slot region and inner belt, with a secondary measurements (6.7 - 35 MeV proton). Such measurements are required to address the following science questions: 1)  Where is the break point in terms of energy of electrons for a given event, below which electrons can be transported into the inner belt from the outer belt but above which electrons cannot, ad what is the injection mechanism? 1)  Where is the break point in terms of energy of electrons for a given event, below which electrons can be transported into the inner belt from the outer belt but above which electrons cannot, ad what is the injection mechanism? 2) What is the CRAND contribution to inner belt electrons, and what is the low energy neutron density near Earth? 3) What is the role of wave-particle interactions in shaping inner-belt electron energy spectra');
        ncwriteatt(Datafilename,"/",'REF1','Mission Overview Paper > Li, X., Kohnert, R., Palo, S., Selesnick, R., Khoo, L.-Y., Schiller, Q., et al. (2022). Two Generations of CubeSat Missions (CSSWE and CIRBE) to Take on the Challenges of Measuring Relativistic Electrons in the Earth’s Magnetosphere. Proceedings of the Small Satellite Conference, SSC22-III-03. https://digitalcommons.usu.edu/smallsat/2022/all2022/152/.')
        ncwriteatt(Datafilename,"/","REF2","Instrument Calibration Paper > Khoo, L.-Y., Li, X., Selesnick, R. S., Schiller, Q., Zhang, K., Zhao, H., et al. (2022). On the challenges of measuring energetic particles in the inner belt: A Geant4 simulation of an energetic particle detector instrument, REPTile-2. Journal of Geophysical Research: Space Physics, 127, e2021JA030249. https://doi.org/10.1029/2021JA030249");
        ncwriteatt(Datafilename,"/","REF3","First Results Paper > Li, X., Selesnick, R., Mei, Y., O'Brien, D., Hogan, B., Xiang, Z., et al. (2024). First results from REPTile-2 measurements onboard CIRBE. Geophysical Research Letters, 51, e2023GL107521. https://doi.org/10.1029/2023GL107521");
        
        % CATDESC
        ncwriteatt(Datafilename,"Epoch","CATDESC","UTC seconds since 2000-01-01 00:00");
        ncwriteatt(Datafilename,"Alt","CATDESC","Altitude");
        ncwriteatt(Datafilename,"Lat","CATDESC","Latitude");
        ncwriteatt(Datafilename,"Lon","CATDESC","Longitude");
        ncwriteatt(Datafilename,"L","CATDESC","L-shell");
        ncwriteatt(Datafilename,"MLT","CATDESC","Magnetic Local Time");
        ncwriteatt(Datafilename,"MLAT","CATDESC","Magnetic latitude");
        ncwriteatt(Datafilename,"Pbins_RNG","CATDESC",'Counts for range proton channels 1:50 (or channels 1:50 for the combined channels)');
        ncwriteatt(Datafilename,"Pbins_PEN","CATDESC",'Counts for penetrating proton channels 1:10 (or channels 51:60 for the combined channels)');
        ncwriteatt(Datafilename,"Ebins_RNG","CATDESC",'Counts for range electron channels 1:50 (or channels 61:110 for the combined channels)');
        ncwriteatt(Datafilename,"Ebins_PEN","CATDESC",'Counts for penetrating electron channels 1:10 (or channels 111:120 for the combined channels)');
        ncwriteatt(Datafilename,"D1","CATDESC",'Counts that trigger D1 (corresponds to channel 121)');
        ncwriteatt(Datafilename,"D2","CATDESC",'Counts that trigger D2 (corresponds to channel 122)');
        ncwriteatt(Datafilename,"D3","CATDESC",'Counts that trigger D3 (corresponds to channel 123)');
        ncwriteatt(Datafilename,"D4","CATDESC",'Counts that trigger D4 (corresponds to channel 124)');
        ncwriteatt(Datafilename,"G","CATDESC",'Counts that trigger guard rings (corresponds to channel 129)');
        ncwriteatt(Datafilename,"D12n","CATDESC",'Counts that trigger D1 but not D2 (corresponds to channel 125)');
        ncwriteatt(Datafilename,"D123n","CATDESC",'Counts that trigger D1, D2, but not D3 (corresponds to channel 126)');
        ncwriteatt(Datafilename,"D1234n","CATDESC",'Counts that trigger D1, D2, D3, but not D4 (corresponds to channel 127)');
        ncwriteatt(Datafilename,"D1234","CATDESC",'Counts that trigger D1, D2, D3, D4 (corresponds to channel 128)');
        ncwriteatt(Datafilename,"IntePrd","CATDESC","Scientific integration period");
        ncwriteatt(Datafilename,"OffAng_timestamp","CATDESC","UTC seconds since 2000-01-01 00:00");
        ncwriteatt(Datafilename,"OffAng","CATDESC","Pointing angle relative to local magnetic field");
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","CATDESC","Flag with values corresponding to pointing mode of spacecraft")
        ncwriteatt(Datafilename,"invalid_data_flag","CATDESC","Flag indicating whether or not data can be considered valid");

        % DISPLAY_TYPE
        ncwriteatt(Datafilename,"Epoch","DISPLAY_TYPE","Time");
        for j = varnames(2:20)
            ncwriteatt(Datafilename,j,"DISPLAY_TYPE","Time series")
        end
        ncwriteatt(Datafilename,"IntePrd","DISPLAY_TYPE","NA");
        ncwriteatt(Datafilename,"OffAng_timestamp","DISPLAY_TYPE","Time");
        ncwriteatt(Datafilename,"OffAng","DISPLAY_TYPE","Time series");
        for j = varnames(24:25)
            ncwriteatt(Datafilename,j,"DISPLAY_TYPE","NA");
        end

        % LABLAXIS
        ncwriteatt(Datafilename,"Epoch","LABLAXIS","Time (UTC seconds since 2000-01-01 00:00)");
        ncwriteatt(Datafilename,"Alt","LABLAXIS","Altitude (km)");
        ncwriteatt(Datafilename,"Lat","LABLAXIS","Latitude (degrees)");
        ncwriteatt(Datafilename,"Lon","LABLAXIS","Longitude (degrees)");
        ncwriteatt(Datafilename,"L","LABLAXIS","L");
        ncwriteatt(Datafilename,"MLT","LABLAXIS","Magnetic Local Time");
        ncwriteatt(Datafilename,"MLAT","LABLAXIS","Magnetic Latitude");
        ncwriteatt(Datafilename,"Pbins_RNG","LABLAXIS","Range Proton Channel Counts");
        ncwriteatt(Datafilename,"Pbins_PEN","LABLAXIS","Penetrating Proton Channel Counts");
        ncwriteatt(Datafilename,"Ebins_RNG","LABLAXIS","Range Electron Channel Counts");
        ncwriteatt(Datafilename,"Ebins_PEN","LABLAXIS","Penetrating Electron Channel Counts");
        for j = varnames(12:20)
            ncwriteatt(Datafilename,j,"LABLAXIS",strjoin([j,"Counts"]))
        end
        ncwriteatt(Datafilename,"IntePrd","LABLAXIS","Integration Period");
        ncwriteatt(Datafilename,"OffAng_timestamp","LABLAXIS","Time (UTC seconds since 2000-01-01 00:00");
        ncwriteatt(Datafilename,"OffAng","LABLAXIS","Pointing Angle (relative to local magnetic field)");
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","LABLAXIS","Spacecraft Pointing Mode")
        ncwriteatt(Datafilename,"invalid_data_flag","LABLAXIS","Data Validity")


        % UNITS
        ncwriteatt(Datafilename,"Epoch","UNITS","UTC Seconds since 2000-01-01 00:00");
        ncwriteatt(Datafilename,"Alt","UNITS","Kilometers");
        ncwriteatt(Datafilename,"Lat","UNITS","Degrees");
        ncwriteatt(Datafilename,"Lon","UNITS","Degrees");
        ncwriteatt(Datafilename,"L","UNITS","None");
        ncwriteatt(Datafilename,"MLT","UNITS","Decimal Hours");
        ncwriteatt(Datafilename,"MLAT","UNITS","Degrees");
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"UNITS","Counts")
        end
        ncwriteatt(Datafilename,"IntePrd","UNITS","Miliseconds");
        ncwriteatt(Datafilename,"OffAng_timestamp","UNITS","Seconds since 2000-01-01 00:00");
        ncwriteatt(Datafilename,"OffAng","UNITS","Degrees");
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","UNITS","0 is nominal science mode, 1 is sun-point static mode, 2 is sun-point star-tracker mode")
        ncwriteatt(Datafilename,"invalid_data_flag","UNITS","0 is valid data, 1 is data that should be used with extreme caution");

        % VALIDMIN
        ncwriteatt(Datafilename,"Epoch","VALIDMIN",single(t(1)));
        ncwriteatt(Datafilename,"Alt","VALIDMIN",0);
        ncwriteatt(Datafilename,"Lat","VALIDMIN",90);
        ncwriteatt(Datafilename,"Lon","VALIDMIN",180);
        ncwriteatt(Datafilename,"L","VALIDMIN",0);
        ncwriteatt(Datafilename,"MLT","VALIDMIN",0);
        ncwriteatt(Datafilename,"MLAT","VALIDMIN",-90);
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"VALIDMIN",0)
        end
        ncwriteatt(Datafilename,"IntePrd","VALIDMIN",999);
        ncwriteatt(Datafilename,"OffAng_timestamp","VALIDMIN",single(thk(1)));
        ncwriteatt(Datafilename,"OffAng","VALIDMIN",0);
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","VALIDMIN",0)
        ncwriteatt(Datafilename,"invalid_data_flag","VALIDMIN",0);

        % VALIDMAX 
        ncwriteatt(Datafilename,"Epoch","VALIDMAX",single(t(num_t)));
        ncwriteatt(Datafilename,"Alt","VALIDMAX",600);
        ncwriteatt(Datafilename,"Lat","VALIDMAX",90);
        ncwriteatt(Datafilename,"Lon","VALIDMAX",180);
        ncwriteatt(Datafilename,"L","VALIDMAX",1.0E3); % above ~12 these numbers have no real physical meaning
        ncwriteatt(Datafilename,"MLT","VALIDMAX",24);
        ncwriteatt(Datafilename,"MLAT","VALIDMAX",90);
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"VALIDMAX",1.25E7) % reached via: 5 second maximum integration period divided by 400 ns minimum time resolution
        end
        ncwriteatt(Datafilename,"IntePrd","VALIDMAX",5.0E3);
        ncwriteatt(Datafilename,"OffAng_timestamp","VALIDMAX",single(thk(num_tHK)));
        ncwriteatt(Datafilename,"OffAng","VALIDMAX",180);
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","VALIDMAX",2)
        ncwriteatt(Datafilename,"invalid_data_flag","VALIDMAX",1);

        % SCALETYP (not in CIRBE combscience but is in minxss)
        for j = varnames(1:7)
            ncwriteatt(Datafilename,j,"SCALETYP","linear")
        end
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"SCALETYP","log")
        end
        for j = varnames(21:25)
            ncwriteatt(Datafilename,j,"SCALETYP","linear")
        end

        % SCALEMIN
        ncwriteatt(Datafilename,"Alt","SCALEMIN",0)
        ncwriteatt(Datafilename,"Lat","SCALEMIN",-90)
        ncwriteatt(Datafilename,"Lon","SCALEMIN",-180)
        ncwriteatt(Datafilename,"L","SCALEMIN",0)
        ncwriteatt(Datafilename,"MLT","SCALEMIN",0)
        ncwriteatt(Datafilename,"MLAT","SCALEMIN",-90)
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"SCALEMIN",0)
        end
        ncwriteatt(Datafilename,"IntePrd","SCALEMIN",999)
        ncwriteatt(Datafilename,"OffAng","SCALEMIN",0)
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","SCALEMIN",-1)
        ncwriteatt(Datafilename,"invalid_data_flag","SCALEMIN",-1);

        % SCALEMAX
        ncwriteatt(Datafilename,"Alt","SCALEMAX",600)
        ncwriteatt(Datafilename,"Lat","SCALEMAX",90)
        ncwriteatt(Datafilename,"Lon","SCALEMAX",180)
        ncwriteatt(Datafilename,"L","SCALEMAX",12) % above ~12 L is physically meaningless
        ncwriteatt(Datafilename,"MLT","SCALEMAX",24)
        ncwriteatt(Datafilename,"MLAT","SCALEMAX",90)
        for j = varnames(8:20)
            ncwriteatt(Datafilename,j,"SCALEMAX",1.0E6)
        end
        ncwriteatt(Datafilename,"IntePrd","SCALEMAX",5.0E3)
        ncwriteatt(Datafilename,"OffAng","SCALEMAX",180)
        ncwriteatt(Datafilename,"spacecraft_pointing_mode","SCALEMAX",3)
        ncwriteatt(Datafilename,"invalid_data_flag","SCALEMAX",2);
    end
end