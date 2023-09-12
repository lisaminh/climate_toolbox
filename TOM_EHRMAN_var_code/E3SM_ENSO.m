function [] = E3SM_ENSO(root,years,Dates,output_dir)

MM = years*12;
ENSO = zeros(MM,1);
time = zeros(MM,2);

DASH = strfind(root,'/');
DD = length(DASH);
MPAS = root(1:DASH(DD-3));
MPAS_file_base = MPAS(DASH(DD-4)+1:DASH(DD-3)-1);

Ocn_Var = 'timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature';

ii = 1;
for yy = min(Dates(:,1)):max(Dates(:,1))
    if yy == min(Dates(:,1))
        m_start = Dates(1,2);
    else
        m_start = 1;
    end
    for mm = m_start:12
        MPAS_file = strcat(MPAS,'/archive/ocn/hist/',MPAS_file_base,'.mpaso.hist.am.timeSeriesStatsMonthly.',num2str(yy),'-',num2str(mm,'%02d'),'-01.nc');
        ENSO(ii) = ncread(MPAS_file,Ocn_Var,[6 1],[1 Inf]);
        time(ii,:) = [yy, mm];
        ii = ii+1;
    end
end

file_out = strcat(output_dir,'/ENSO.nc');

nccreate(file_out,'ENSO','Dimensions',{'time',MM});
nccreate(file_out,'time','Dimensions',{'time',MM,'date',2});

ncwrite(file_out,'ENSO',ENSO);
ncwrite(file_out,'time',time);