clear
clc

% Root directory for E3SM output
root='/pscratch/sd/l/lmnguyen/TEM_data/'

%root='/pscratch/sd/w/wagmanbe/E3SM_simulations/CLDERA/v2.LR.WCYCL20TR.pmcpu.limvar.ens1/archive/atm/hist/';
%root='C:/Users/lisa/Downloads/variability-master-p1_lim_v_full-MATLAB/variability-master-p1_lim_v_full-MATLAB/p1_lim_v_full/MATLAB/';
% Name of first E3SM output file
first_file='v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES.nc'

% Number of complete years simulated
years=1;

% Directory for output files
output_dir='/global/homes/l/lmnguyen/CLDERA/data/TEM_data';
%output_dir='/pscratch/sd/t/tsehrma/E3SM_QOI_ens';
%output_dir='C:/Users/lisa/Downloads/variability-master-p1_lim_v_full-MATLAB/variability-master-p1_lim_v_full-MATLAB/p1_lim_v_full/MATLAB';

% Use file directory (stores transformation matrices within file)
%Use_files = '/global/homes/t/tsehrma/E3SM_QOI/Use_files';
Use_files='/global/homes/l/lmnguyen/CLDERA/data/Sph_ZM_matrices';

% used in E3SM_find_dates to put the file_info.mat should really put this in that script
if exist(strcat(output_dir,'/run_files'),'dir')==0
    mkdir(strcat(output_dir,'/run_files'))
end

[Dates,file_base] = E3SM_find_dates(root,first_file,years,output_dir);

file_check = strcat(Use_files,'/ZM.nc');
%TODO: check for specific ZM_matrix file by renaming matrix files in a certain convention
if exist(file_check,'file')==0
    lat = read_E3SM_QOI(file_base,'1',Dates(1,:),'lat');
    lat_ZM = (-89.5:.5:89.5)';
    L = 150;
    
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

%TODO: change this back
Dates_full=Dates;
%Dates=Dates_full(1:3,:);
disp('Calculating TEM variables')
[PSI, divF, Flat,Fp,V_TEM, OMEGA_TEM, W_TEM, lat,lev, time] = E3SM_TEM(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir);


%file_out = strcat(output_dir,'/TEMv4_vars.nc');
%ncid = netcdf.create(file_out, 'NOCLOBBER')
%dimid(1)=netcdf.defDim(ncid,'time',TTT)
%dimid(2)=netcdf.defDim(ncid,'lev',LL)
%dimid(3)=netcdf.defDim(ncid,'lat',MM)
%tim_ID=netcdf.defVar(ncid,'time','double',dimidtim);
%
%varid=netcdf.defVar(ncid,'F_lat', '', dimid
%
%
%%
%nccreate(file_out,'F_lat','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
%nccreate(file_out,'F_p','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
%nccreate(file_out,'OMEGA_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
%nccreate(file_out,'V_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
%nccreate(file_out,'W_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
%%
%ncwrite(file_out,'F_lat',Flat);
%ncwrite(file_out,'F_p',Fp);
%ncwrite(file_out,'OMEGA_TEM',OMEGA_TEM);
%ncwrite(file_out,'V_TEM',V_TEM);
%ncwrite(file_out,'W_TEM',W_TEM);
