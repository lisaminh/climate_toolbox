function [Dates,file_base] = E3SM_find_dates(root,first_file,years,output_dir)
% This function returns all of the start dates for each file name in an E3SM output 

%clear
%clc
%config

%# Root directory for E3SM output
%root='/global/cscratch1/sd/wagmanbe/E3SM_simulations/CLDERA/v2.LR.WCYCL20TR_0161_ne30pg2_prog_cori/archive/atm/hist/';
%
%# Name of first output file
%first_file='v2.LR.WCYCL20TR_0161_ne30pg2_prog_cori.eam.h1.1850-01-01-00000.nc';
%
%# Number of complete years simulated
%years=165;

%# Directory for output files
%output_dir='/global/cscratch1/sd/tsehrma/QOI_test';

DOTS = strfind(first_file,'.');
DD = length(DOTS);

file_base = strcat(root,first_file(1:DOTS(DD-1)-2));

year1 = str2double(first_file(DOTS(DD-1)+1:DOTS(DD-1)+4));
month = str2double(first_file(DOTS(DD-1)+6:DOTS(DD-1)+7));
day = str2double(first_file(DOTS(DD-1)+9:DOTS(DD-1)+10));
% month = 1;
% day = 1;

FF = ceil(years*365/30);

month_day = [31 28 31 30 31 30 31 31 30 31 30 31];

Dates = [year1*ones(FF,1),ones(FF,1),(0:FF-1)'*30+1+sum(month_day(1:month-1))];
for dd = 1:FF
    while Dates(dd,3) > 365
        Dates(dd,3) = Dates(dd,3) - 365;
        Dates(dd,1) = Dates(dd,1) + 1;
    end
    while Dates(dd,3) > month_day(Dates(dd,2))
        Dates(dd,3) = Dates(dd,3) - month_day(Dates(dd,2));
        Dates(dd,2) = Dates(dd,2) + 1;
    end
end

out_file = strcat(output_dir,'/run_files/file_info.mat');
save(out_file,'Dates','file_base')
