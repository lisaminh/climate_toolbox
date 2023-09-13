% Calculates EP flux components and flux divergence from E3SM output

%function [F_theta,F_p,divF,Fp_trop,divF_strat] = E3SM_EPflux(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir)
function [Flat,Fp, V_TEM, OMEGA_TEM, W_TEM, lat_ZM,lev,time] = E3SM_TEM(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir)

h_num = '1';

kappa = 0.287;
Re = 6.3781e6;
Om = 7.2921e-5;

%DD = length(Dates, dim=1);
DD = size(Dates, 1);

%TODO: change this back! assumes there are 30 days in each file
TTT=30*DD;
time=1:TTT;
%TTT = 365*years; % number of days
%TTT = floor(365*years/7)+1;
H=7000; %scale height set to 7km
dlat = (lat_ZM(2)-lat_ZM(1))*pi/180;
MM = length(lat_ZM); % #lats
LL = 72; %TODO: manually set levels, have this be an input parameter
Fp = zeros(MM,LL,TTT);
Fp = zeros(MM,LL,TTT);
divF = zeros(MM,LL,TTT);
V_TEM = zeros(MM,LL,TTT);
OMEGA_TEM = zeros(MM,LL,TTT);
W_TEM = zeros(MM,LL,TTT);

for dd = 1:DD
    U = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'U');
    V = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'V');
    T = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'T');
    Omega = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'OMEGA');

    if dd == 1
        lev = read_E3SM_QOI(file_base,h_num,Dates(1,:),'lev');
        [NN LL TT] = size(T); % #ncol, #lev, #time
                
        f = 2*Om*repmat(sind(lat_ZM),[1 LL TT]);
        coslat = repmat(cosd(lat_ZM),[1 LL TT]);
    elseif dd == DD
	    TT_old = TT;
	    [~,~,TT] = size(T);
        f = 2*Om*repmat(sind(lat_ZM),[1 LL TT]);
        coslat = repmat(cosd(lat_ZM),[1 LL TT]);
    end
    
    Theta = (1000^kappa)*T.*repmat(lev'.^-kappa,[NN 1 TT]); %here using lev, not pressure (dcan take top 33 levs where hybm=0, or up to lev=170.744079521582)
    
    %calculate theta_bar in the new lats
    Theta_b = reshape(ZM*reshape(Theta,[NN LL*TT]),[MM LL TT]);
    % calculate dtheta_bar/dp using centered difference
    dTheta_bdp = [(Theta_b(:,2,:)-Theta_b(:,1,:))/(lev(2)-lev(1)), (Theta_b(:,3:LL,:)-Theta_b(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (Theta_b(:,LL,:)-Theta_b(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    
    %compute zonal mean in the natural coordinates
    U_bar_nat = reshape(ZM_nat*reshape(U,[NN LL*TT]),[NN LL TT]);
    V_bar_nat = reshape(ZM_nat*reshape(V,[NN LL*TT]),[NN LL TT]);
    Theta_bar_nat = reshape(ZM_nat*reshape(Theta,[NN LL*TT]),[NN LL TT]);
    Omega_bar_nat = reshape(ZM_nat*reshape(Omega,[NN LL*TT]),[NN LL TT]);
    
    %take the eddy components
    U_p = U - U_bar_nat;
    V_p = V - V_bar_nat;
    Theta_p = Theta - Theta_bar_nat;
    Omega_p = Omega - Omega_bar_nat;
    
    %take the quadratic eddy components in the new lats
    VpThetap = reshape(ZM*reshape(V_p.*Theta_p,[NN LL*TT]),[MM LL TT]);
    UpOmegap = reshape(ZM*reshape(U_p.*Omega_p,[NN LL*TT]),[MM LL TT]);
    UpVp = reshape(ZM*reshape(U_p.*V_p,[NN LL*TT]),[MM LL TT]);
    %VpTp = reshape(ZM*reshape(V_p.*Theta_p,[NN LL*TT]),[MM LL TT])./dTheta_bdp;
    
    %F_theta = -Re*coslat.*UpVp;
    %F_p = Re*coslat.*f.*VpTp;
    
    %Calculate eddy stream function Psi
    Psi = VpThetap./dTheta_bdp;
    dPsidp=[(Psi(:,2,:)-Psi(:,1,:))/(lev(2)-lev(1)), (Psi(:,3:LL,:)-Psi(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (Psi(:,LL,:)-Psi(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    Psi_cos = Psi.*coslat;
    dPsi_cosdlat = [(Psi_cos(1,:,:)-Psi_cos(2,:,:))/dlat; (Psi_cos(1:MM-2,:,:)-Psi_cos(3:MM,:,:))/(2*dlat); (Psi_cos(MM-1,:,:)-Psi_cos(MM,:,:))/dlat];

    %calculateTEM northward and vertical velocities
    Omega_b = reshape(ZM*reshape(Omega,[NN LL*TT]),[MM LL TT]); %Omega_bar in new lats
    V_b = reshape(ZM*reshape(V,[NN LL*TT]),[MM LL TT]); %V_bar in new lats

    Vtem=V_b-dPsidp; %in m/s
    Omegatem = Omega_b+(dPsi_cosdlat./coslat)/Re; %in hPa/s
    Wtem = (Omegatem./repmat(lev',[MM,1,TT]))*H/100; %in m/s

    %calculate Eliassen-Palm flux
    U_b = reshape(ZM*reshape(U,[NN LL*TT]),[MM LL TT]); %U_bar in new lats
    dUdp = [(U_b(:,2,:)-U_b(:,1,:))/(lev(2)-lev(1)), (U_b(:,3:LL,:)-U_b(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (U_b(:,LL,:)-U_b(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    
    U_b_cos = U_b.*coslat;
    dU_b_cosdlat = [(U_b_cos(1,:,:)-U_b_cos(2,:,:))/dlat; (U_b_cos(1:MM-2,:,:)-U_b_cos(3:MM,:,:))/(2*dlat); (U_b_cos(MM-1,:,:)-U_b_cos(MM,:,:))/dlat];
    
    %Original ones assumed Psi=VpTp
    %F_lat = Re*coslat.*(-1*UpVp+dUdp.*VpTp);
    %F_p = Re*coslat.*VpTp.*(f-((Re*coslat).^-1).*dU_b_cosdlat);

    F_lat = Re*coslat.*(-1*UpVp+dUdp.*Psi);
    F_p = Re*coslat.* ( (f-((Re*coslat).^-1).*dU_b_cosdlat).*Psi-UpOmegap);
    
    F_theta_coslat = F_lat.*coslat;
    dFdlat = ((Re*coslat).^-1).*[(F_theta_coslat(1,:,:)-F_theta_coslat(2,:,:))/dlat; (F_theta_coslat(1:MM-2,:,:)-F_theta_coslat(3:MM,:,:))/(2*dlat); (F_theta_coslat(MM-1,:,:)-F_theta_coslat(MM,:,:))/dlat];
    dFdp = [(F_p(:,2,:)-F_p(:,1,:))/(lev(2)-lev(1)), (F_p(:,3:LL,:)-F_p(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (F_p(:,LL,:)-F_p(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    %divF = (dFdlat)./(Re*coslat) + dFdp;
    
    if dd == DD
   	    TT_use = TTT-((dd-1)*TT_old);
        Fp(:,:,(dd-1)*TT_old+1:TTT)=F_p;
        Flat(:,:,(dd-1)*TT_old+1:TTT)=F_lat;
        V_TEM(:,:,(dd-1)*TT_old+1:TTT)=Vtem;
        OMEGA_TEM(:,:,(dd-1)*TT_old+1:TTT)=Omegatem;
        W_TEM(:,:,(dd-1)*TT_old+1:TTT)=Wtem;

    else
        Fp(:,:,(dd-1)*TT+1:dd*TT) = F_p;
        Flat(:,:,(dd-1)*TT+1:dd*TT) = F_lat;
        V_TEM(:,:,(dd-1)*TT+1:dd*TT) = Vtem;
        OMEGA_TEM(:,:,(dd-1)*TT+1:dd*TT) = Omegatem;
        W_TEM(:,:,(dd-1)*TT+1:dd*TT) = Wtem;
    end

end
file_out = strcat(output_dir,'/TEM_vars.nc');

nccreate(file_out,'F_lat','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
nccreate(file_out,'F_p','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
nccreate(file_out,'OMEGA_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
nccreate(file_out,'V_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});
nccreate(file_out,'W_TEM','Dimensions',{'lat', MM,'lev',LL ,'time',TTT});

% For some reason when I include the coordinates, they do not add
%nccreate(file_out, 'lat', 'Dimensions',{'lat', MM});
%nccreate(file_out, 'time', 'Dimensions',{'time', TTT});
%nccreate(file_out, 'lev', 'Dimensions',{'lev', LL});

ncwrite(file_out,'F_lat',Flat);
ncwrite(file_out,'F_p',Fp);
ncwrite(file_out,'OMEGA_TEM',OMEGA_TEM);
ncwrite(file_out,'V_TEM',V_TEM);
ncwrite(file_out,'W_TEM',W_TEM);
%ncwrite(file_out,'lat',lat_ZM);
%ncwrite(file_out,'time',time);
%ncwrite(file_out,'lev',lev);

