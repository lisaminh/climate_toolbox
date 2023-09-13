function [SAMi,time] = E3SM_SAMi(SAM,Dates,file_base,years,output_dir)
% This function will calculate the NAM index for a given NAM base`

h_num = '1';
FF = length(Dates);

time = Dates(1,1)+(0:years*365);
TT = length(time);
% TT = (max(Dates(:,1))-min(Dates(:,1))+1)*365;
%TT = floor((max(Dates(:,1))-min(Dates(:,1))+1)*365/7)+1;

ff = 1;
lat = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'lat');
lon = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'lon');
time = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'time');
TTT = length(time);

Use = lat<-20;
NN = sum(Use);

Z = zeros(NN,TT);
Z_temp = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z700');
%Z_temp = squeeze(read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z3',[1 52 1],[Inf 1 Inf]));
Z(:,1:TTT) = Z_temp(Use,:);

for ff = 2:FF-1
    Z_temp = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z700');
%    Z_temp = squeeze(read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z3',[1 52 1],[Inf 1 Inf]));
    Z(:,(ff-1)*TTT+1:ff*TTT) = Z_temp(Use,:);
end
Z_temp = read_E3SM_QOI(file_base,h_num,Dates(FF,:),'Z700');
%Z_temp = squeeze(read_E3SM_QOI(file_base,h_num,Dates(FF,:),'Z3',[1 52 1],[Inf 1 Inf]));
Z(:,(FF-1)*TTT+1:TT) = Z_temp(Use,1:TT-(FF-1)*TTT);

SAMi = ((Z-repmat(mean(Z,2),[1 TT]))./repmat(std(Z,1,2),[1 TT]))'*SAM;

file_out = strcat(output_dir,'/SAM_index.nc');

time = (0:TT-1);

nccreate(file_out,'SAMi','Dimensions',{'time',TT})
nccreate(file_out,'time','Dimensions',{'time',TT})

ncwrite(file_out,'SAMi',SAMi);
ncwrite(file_out,'time',time);
