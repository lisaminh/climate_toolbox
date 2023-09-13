clear
clc

% Root directory for E3SM output
root='/pscratch/sd/w/wagmanbe/E3SM_simulations/CLDERA/v2.LR.WCYCL20TR.0211.trc.pmcpu/archive/atm/hist/';

% Name of first E3SM output file
first_file='v2.LR.WCYCL20TR.0211.trc.pmcpu.eam.h1.1850-01-01-00000.nc';

% Number of complete years simulated
years=165;

% Directory for output files
output_dir='/pscratch/sd/t/tsehrma/E3SM_QOI_ens';

% Use file directory (stores transformation matrices within file)
Use_files = '/global/homes/t/tsehrma/E3SM_QOI/Use_files';

if exist(strcat(output_dir,'/run_files'),'dir')==0
    mkdir(strcat(output_dir,'/run_files'))
end

[Dates,file_base] = E3SM_find_dates(root,first_file,years,output_dir);

disp('Calculating ENSO')
E3SM_ENSO(root,years,Dates,output_dir);

disp('Calculating QBO')
E3SM_QBO(Dates,file_base,years,output_dir);

lev = 25;
E3SM_QBO_lev(Dates,file_base,years,output_dir,lev);

disp('Calculating Polar Vortex')
E3SM_Polar_Vortex(Dates,file_base,years,output_dir);

if exist(strcat(Use_files,'/NAM.nc'),'file')==0
    year1 = 1950;
    [NAM,lat_use,lon_use,Use] = E3SM_NAM_build(Dates,file_base,years,Use_files,year1);
else
    NAM = ncread(strcat(Use_files,'/NAM.nc'),'NAM');
end
disp('Calculating NAM')
[NAMi, time] = E3SM_NAMi(NAM,Dates,file_base,years,output_dir);

if exist(strcat(Use_files,'/SAM.nc'),'file')==0
    year1 = 1950;
    [SAM,lat_use2,lon_use2,Use2] = E3SM_SAM_build(Dates,file_base,years,Use_files,year1);
else
    SAM = ncread(strcat(Use_files,'/SAM.nc'),'SAM');
end
disp('Calculating SAM')
[SAMi, time] = E3SM_SAMi(SAM,Dates,file_base,years,output_dir);

file_check = strcat(Use_files,'/ZM.nc');
if exist(file_check,'file')==0
    lat = read_E3SM_QOI(file_base,'1',Dates(1,:),'lat');
    lat_ZM = (-89.5:.5:89.5);
    L = 100;
    
    [ZM, ZM_nat] = Sph_Zonal_Mean(lat,lat_ZM,L);

    nccreate(file_check,'ZM','Dimensions',{'lat_ZM',length(lat_ZM),'lat',length(lat)});
    nccreate(file_check,'lat_ZM','Dimensions',{'lat_ZM',length(lat_ZM)});
    nccreate(file_check,'ZM_nat','Dimensions',{'lat',length(lat),'lat',length(lat)});
    nccreate(file_check,'lat','Dimensions',{'lat',length(lat)});

    ncwrite(file_check,'ZM',ZM);
    ncwrite(file_check,'lat_ZM',lat_ZM);
    ncwrite(file_check,'ZM_nat',ZM_nat);
    ncwrite(file_check,'lat',lat);
else
    ZM = ncread(file_check,'ZM');
    lat_ZM = ncread(file_check,'lat_ZM');
    ZM_nat = ncread(file_check,'ZM_nat');
    lat = ncread(file_check,'lat');
end
 
disp('Calculating EPflux')
[Fp_surf,Fp_trop,divF_strat] = E3SM_EPflux(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir);

exit
