function [NAM,lat_use,lon_use,Use] = E3SM_NAM_build(Dates,file_base,years,output_dir,year1)
% This is the function builds the NAM index

h_num = '1';
FF = length(Dates);

TT = (max(Dates(:,1))-year1+1)*365;

f_start = find(Dates(:,1)==year1-1,1,'last');
lat = read_E3SM_QOI(file_base,h_num,Dates(f_start,:),'lat');
lon = read_E3SM_QOI(file_base,h_num,Dates(f_start,:),'lon');
time = read_E3SM_QOI(file_base,h_num,Dates(f_start,:),'time');
TTT = length(time);

Use = lat>20;
NN = sum(Use);

time = time-floor(time/365)*365;
t_start = find(time==0);
t1 = TTT-t_start+1;

Z = zeros(NN,TT);
Z_temp = read_E3SM_QOI(file_base,h_num,Dates(f_start,:),'Z1000',[1 t_start],[Inf Inf]);
Z(:,1:t1) = Z_temp(Use,:);

for ff = f_start+1:FF-1
    Z_temp = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z1000');
    Z(:,(ff-f_start-1)*TTT+t1+1:(ff-f_start)*TTT+t1) = Z_temp(Use,:);
end
Z_temp = read_E3SM_QOI(file_base,h_num,Dates(FF,:),'Z1000');
Z(:,(FF-f_start-1)*TTT+t1+1:TT) = Z_temp(Use,1:TT-((FF-f_start-1)*TTT+t1));

A = 91;
B = (1/A)*ones(1,A);
Z_use = [Z(:,TT-44:TT), Z, Z(:,1:45)];
Z_use = conv2(1,B,Z_use,'valid');

time = (0:TT-1);
time = time-floor(time/365)*365;
T_use = or(time<=120,time>304);

Z_use = Z_use(:,T_use);

Z_use = Z_use - repmat(mean(Z_use,2),[1,sum(T_use)]);

covZ = (1/sum(T_use))*Z_use*Z_use';

[EOF,Val] = eig(covZ,'vector');
[M,I] = max(abs(Val));

lat_use = lat(Use);
lon_use = lon(Use);

%if I == 1
if EOF(find(lat_use>85,1,'first'),1)>0
    NAM = -1*EOF(:,I);
else
    NAM = EOF(:,I);
end

file_out = strcat(output_dir,'/NAM.nc');

nccreate(file_out,'NAM','Dimensions',{'N_use',NN});
nccreate(file_out,'lat_use','Dimensions',{'N_use',NN});
nccreate(file_out,'lon_use','Dimensions',{'N_use',NN});
nccreate(file_out,'N_use','Dimensions',{'N',length(Use)});

ncwrite(file_out,'NAM',NAM);
ncwrite(file_out,'lat_use',lat_use);
ncwrite(file_out,'lon_use',lon_use);
ncwrite(file_out,'N_use',double(Use));

%NAM_i = ((Z-repmat(mean(Z,2),[1 TT]))./repmat(std(Z,1,2),[1 TT]))'*NAM;
%
%figure
%plot((0:TT-1),NAM_i)
%
%R = 6373;
%X = R*cosd(lat(Use)).*sind(lon(Use));
%Y = -R*cosd(lat(Use)).*cosd(lon(Use));
%
%cmap = [0 0 1;
%        0 .33 1;
%        0 .66 1;
%        0 1 1;
%        0 1 0;
%        1 1 0;
%        1 .66 0;
%        1 .33 0;
%        1 0 0;
%        .5 0 0];
%
%CC = length(cmap);
%
%c_max = max(NAM);
%c_min = min(NAM);
%c_int = (c_max-c_min)/CC;
%cons = (c_min:c_int:c_max);
%
%figure
%hold on
%for cc = 1:CC
%    Data_use = and(NAM>=cons(cc),NAM<cons(cc+1));
%    plot(X(Data_use),Y(Data_use),'.','Markersize',4,'color',cmap(cc,:))
%end
%hold off
%box on
%colormap(cmap)
%caxis([c_min c_max])
%colorbar
