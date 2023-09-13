function [] = E3SM_Polar_Vortex(Dates,file_base,years,output_dir)
% Derives and saves the SSW and Polar Vortex Moment time series for a given E3SM run

FF = length(Dates);

h_num = '1';

ff = 1;
lat = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'lat');
lon = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'lon');
area = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'area');
U = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'U010');
Z = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z010');
[~, TT] = size(U);

Use_N = and(lat<65,lat>60);
Use_S = and(lat<-60,lat>-65);

Weight_N = repmat(area(Use_N)/sum(area(Use_N)),[1 TT]);
Weight_S = repmat(area(Use_S)/sum(area(Use_S)),[1 TT]);

time = Dates(1,1)+(0:years*365);
%time = Dates(1,1)+(0:7:years*365);
TTT = length(time);

SSW = zeros(TTT,2);
Z_edge = zeros(TTT,2);

SSW(1:TT,1) = squeeze(sum(U(Use_N,:,:).*Weight_N,1));
SSW(1:TT,2) = squeeze(sum(U(Use_S,:,:).*Weight_S,1));
Z_edge(1:TT,1) = squeeze(sum(Z(Use_N,:).*Weight_N,1));
Z_edge(1:TT,2) = squeeze(sum(Z(Use_S,:).*Weight_S,1));
for ff = 2:FF-1
    U = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'U010');
    SSW((ff-1)*TT+1:TT*ff,1) = squeeze(sum(U(Use_N,:).*Weight_N,1));
    SSW((ff-1)*TT+1:TT*ff,2) = squeeze(sum(U(Use_S,:).*Weight_S,1));
    
    Z = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z010');
    Z_edge((ff-1)*TT+1:TT*ff,1) = squeeze(sum(Z(Use_N,:).*Weight_N,1));
    Z_edge((ff-1)*TT+1:TT*ff,2) = squeeze(sum(Z(Use_S,:).*Weight_S,1));
end
U = read_E3SM_QOI(file_base,h_num,Dates(FF,:),'U010');
SSW((FF-1)*TT+1:TTT,1) = squeeze(sum(U(Use_N,1:TTT-(FF-1)*TT).*Weight_N(:,1:TTT-(FF-1)*TT),1));
SSW((FF-1)*TT+1:TTT,2) = squeeze(sum(U(Use_S,1:TTT-(FF-1)*TT).*Weight_S(:,1:TTT-(FF-1)*TT),1));

Z = read_E3SM_QOI(file_base,h_num,Dates(FF,:),'Z010');
Z_edge((FF-1)*TT+1:TTT,1) = squeeze(sum(Z(Use_N,1:TTT-(FF-1)*TT).*Weight_N(:,1:TTT-(FF-1)*TT),1));
Z_edge((FF-1)*TT+1:TTT,2) = squeeze(sum(Z(Use_S,1:TTT-(FF-1)*TT).*Weight_S(:,1:TTT-(FF-1)*TT),1));

file_out = strcat(output_dir,'/Polar_Vortex.nc');

nccreate(file_out,'SSW','Dimensions',{'time',TTT,'Hemisphere',2});
nccreate(file_out,'Z_edge','Dimensions',{'time',TTT,'Hemisphere',2});
nccreate(file_out,'time','Dimensions',{'time',TTT});

ncwrite(file_out,'SSW',SSW);
ncwrite(file_out,'Z_edge',Z_edge);
ncwrite(file_out,'time',time);

%-----------------

Re = 6.357e3;

%A = 4*30;
A = 4;
B = (1/A)*ones(A,1);

Z_NH_edge = Z_edge(:,1)/1e3;
Z_NH_edge = [Z_NH_edge(TT-A/2+1:TT); Z_NH_edge; Z_NH_edge(1:A/2-1)];
Z_NH_edge = conv(Z_NH_edge,B,'valid');

Z_SH_edge = Z_edge(:,2)/1e3;
Z_SH_edge = [Z_SH_edge(TT-A/2+1:TT); Z_SH_edge; Z_SH_edge(1:A/2-1)];
Z_SH_edge = conv(Z_SH_edge,B,'valid');

AR = zeros(TTT,2);
Xc = zeros(TTT,2);
Yc = zeros(TTT,2);

area = (4/3)*pi*(Re^3)/sum(area)*area;
TT_use = TT;
     
NH = lat>0;
SH = lat<0;

X_NH = Re*sqrt(2*(1-sind(lat(NH)))).*cosd(lon(NH));
Y_NH = Re*sqrt(2*(1-sind(lat(NH)))).*sind(lon(NH));
X_SH = Re*sqrt(2*(1-sind(lat(SH)))).*cosd(lon(SH));
Y_SH = Re*sqrt(2*(1-sind(lat(SH)))).*sind(lon(SH));

for ff = 1:FF
    Z = read_E3SM_QOI(file_base,h_num,Dates(ff,:),'Z010');
    if ff == FF
        [~, TT_use] = size(Z);
    end
    
    Z = Z/1e3; % convert to km
    Z_NH = Z(NH,:);
    Z_SH = Z(SH,:);
    
    for tt = 1:TT_use
        IN = Z_NH(:,tt)<Z_NH_edge((ff-1)*TT+tt);
        m00 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*area(IN));
        m10 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*X_NH(IN).*area(IN));
        m01 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*Y_NH(IN).*area(IN));
        Xc((ff-1)*TT+tt,1) = m10/m00;
        Yc((ff-1)*TT+tt,1) = m01/m00;
        
        J20 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*((X_NH(IN)-Xc((ff-1)*TT+tt,1)).^2).*area(IN));
        J02 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*((Y_NH(IN)-Yc((ff-1)*TT+tt,1)).^2).*area(IN));
        J11 = sum((Z_NH(IN,tt)-Z_NH_edge((ff-1)*TT+tt)).*(X_NH(IN)-Xc((ff-1)*TT+tt,1)).*(Y_NH(IN)-Yc((ff-1)*TT+tt,1)).*area(IN));
        AR((ff-1)*TT+tt,1) = sqrt(abs((J20+J02+sqrt(4*J11^2+(J20-J02)^2))/(J20+J02-sqrt(4*J11^2+(J20-J02)^2))));
        
        IN = Z_SH(:,tt)<Z_SH_edge((ff-1)*TT+tt);
        m00 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*area(IN));
        m10 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*X_SH(IN).*area(IN));
        m01 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*Y_SH(IN).*area(IN));
        Xc((ff-1)*TT+tt,2) = m10/m00;
        Yc((ff-1)*TT+tt,2) = m01/m00;
        
        J20 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*((X_SH(IN)-Xc((ff-1)*TT+tt,2)).^2).*area(IN));
        J02 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*((Y_SH(IN)-Yc((ff-1)*TT+tt,2)).^2).*area(IN));
        J11 = sum((Z_SH(IN,tt)-Z_SH_edge((ff-1)*TT+tt)).*(X_SH(IN)-Xc((ff-1)*TT+tt,2)).*(Y_SH(IN)-Yc((ff-1)*TT+tt,2)).*area(IN));
        AR((ff-1)*TT+tt,2) = sqrt(abs((J20+J02+sqrt(4*J11^2+(J20-J02)^2))/(J20+J02-sqrt(4*J11^2+(J20-J02)^2))));
    end
end

Latc = asind(sqrt(1-((Xc.^2 + Yc.^2)/(2*Re^2))));
Lonc = atand(Yc./Xc);
Lonc(Xc<0) = Lonc(Xc<0)+180;
Lonc(Lonc<0) = Lonc(Lonc<0)+360;
AR(AR<.2) = NaN;
AR(AR<1) = AR(AR<1).^-1;

nccreate(file_out,'Latc','Dimensions',{'time',TTT,'Hemisphere',2});
nccreate(file_out,'Lonc','Dimensions',{'time',TTT,'Hemisphere',2});
nccreate(file_out,'AR','Dimensions',{'time',TTT,'Hemisphere',2});

ncwrite(file_out,'Latc',Latc);
ncwrite(file_out,'Lonc',Lonc);
ncwrite(file_out,'AR',AR);
