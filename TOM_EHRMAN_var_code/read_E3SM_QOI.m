function Data = read_E3SM_QOI(file_base,h_num,Date,Var,start,count)
% reads E3SM data

file_tail = '-00000.nc';

file_name = strcat(file_base,h_num,'.',num2str(Date(1),'%04d'),'-',num2str(Date(2),'%02d'),'-',num2str(Date(3),'%02d'),file_tail);

if nargin > 4
    Data = ncread(file_name,Var,start,count);
else
    Data = ncread(file_name,Var);
end
