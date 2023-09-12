function [] = E3SM_QBO(Dates,file_base,years,output_dir)
% Derives and saves the QBO time series for a given E3SM run

FF = length(Dates);

h_num = '1';

ff = 1;
lat = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'lat');
area = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'area');
U = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'U050');
[~, TT] = size(U);

Use = abs(lat) <= 5;
Weight = repmat(area(Use)/sum(area(Use)),[1 TT]);

time = Dates(1,1)+(0:years*365);
%time = Dates(1,1)+(0:7:years*365);
TTT = length(time);

QBO = zeros(TTT,1);

QBO(1:TT) = squeeze(sum(U(Use,:,:).*Weight,1));
for ff = 2:FF-1
    U = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'U050');
    QBO((ff-1)*TT+1:TT*ff) = squeeze(sum(U(Use,:).*Weight,1));
end
U = read_E3SM_QOI(file_base,h_num,Dates(FF,:),'U050');
QBO((FF-1)*TT+1:TTT) = squeeze(sum(U(Use,1:TTT-(FF-1)*TT).*Weight(:,1:TTT-(FF-1)*TT),1));

file_out = strcat(output_dir,'/QBO_50.nc');

nccreate(file_out,'QBO','Dimensions',{'time',TTT});
nccreate(file_out,'time','Dimensions',{'time',TTT});

ncwrite(file_out,'QBO',QBO);
ncwrite(file_out,'time',time);
