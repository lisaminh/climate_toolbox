% Calculates EP flux components and flux divergence from E3SM output

%function [F_theta,F_p,divF,Fp_trop,divF_strat] = E3SM_EPflux(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir)
function [Fp_surf,Fp_trop,divF_strat] = E3SM_EPflux(ZM,ZM_nat,lat_ZM,lat,Dates,file_base,years,output_dir)

h_num = '1';

kappa = 0.287;
Re = 6.3781e6;
Om = 7.2921e-5;

DD = length(Dates);

TTT = 365*years;
%TTT = floor(365*years/7)+1;

Fp_surf = zeros(TTT,1);
Fp_trop = zeros(TTT,1);
divF_strat = zeros(TTT,1);

for dd = 1:DD
    U = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'U');
    V = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'V');
    T = read_E3SM_QOI(file_base,h_num,Dates(dd,:),'T');
    
    if dd == 1
        lev = read_E3SM_QOI(file_base,h_num,Dates(1,:),'lev');
        [NN LL TT] = size(T);
        MM = length(lat_ZM);
        
        M_20 = find(lat_ZM==20);
        M_45 = find(lat_ZM==45);
        M_65 = find(lat_ZM==65);
        
	M_n20 = find(lat_ZM==-20);
        M_n45 = find(lat_ZM==-45);
        M_n65 = find(lat_ZM==-65);
        
	L_850 = 58;
        L_200 = 35;
        L_10 = 13;
        
        dlat = (lat_ZM(2)-lat_ZM(1))*pi/180;
        
        f = 2*Om*repmat(sind(lat_ZM),[1 LL TT]);
        coslat = repmat(cosd(lat_ZM),[1 LL TT]);
    elseif dd == DD
	TT_old = TT;
	[~,~,TT] = size(T);
        f = 2*Om*repmat(sind(lat_ZM),[1 LL TT]);
        coslat = repmat(cosd(lat_ZM),[1 LL TT]);
    end
    
    Theta = (1000^kappa)*T.*repmat(lev'.^-kappa,[NN 1 TT]);
    
    %U = squeeze(mean(U,3));
    %V = squeeze(mean(V,3));
    %Theta = squeeze(mean(Theta,3));
    
    Theta_b = reshape(ZM*reshape(Theta,[NN LL*TT]),[MM LL TT]);
    dTheta = [(Theta_b(:,2,:)-Theta_b(:,1,:))/(lev(2)-lev(1)), (Theta_b(:,3:LL,:)-Theta_b(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (Theta_b(:,LL,:)-Theta_b(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    
    U_p = U - reshape(ZM_nat*reshape(U,[NN LL*TT]),[NN LL TT]);
    V_p = V - reshape(ZM_nat*reshape(V,[NN LL*TT]),[NN LL TT]);
    Theta_p = Theta - reshape(ZM_nat*reshape(Theta,[NN LL*TT]),[NN LL TT]);
    
    UpVp = reshape(ZM*reshape(U_p.*V_p,[NN LL*TT]),[MM LL TT]);
    VpTp = reshape(ZM*reshape(V_p.*Theta_p,[NN LL*TT]),[MM LL TT])./dTheta;
    
    %F_theta = -Re*coslat.*UpVp;
    %F_p = Re*coslat.*f.*VpTp;
    
    U_b = reshape(ZM*reshape(U,[NN LL*TT]),[MM LL TT]);
    dUdp = [(U_b(:,2,:)-U_b(:,1,:))/(lev(2)-lev(1)), (U_b(:,3:LL,:)-U_b(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (U_b(:,LL,:)-U_b(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    
    U_b = U_b.*coslat;
    dUdlat = [(U_b(1,:,:)-U_b(2,:,:))/dlat; (U_b(1:MM-2,:,:)-U_b(3:MM,:,:))/(2*dlat); (U_b(MM-1,:,:)-U_b(MM,:,:))/dlat];
    
    F_theta = Re*coslat.*(-1*UpVp+dUdp.*VpTp);
    F_p = Re*coslat.*VpTp.*(f-((Re*coslat).^-1).*dUdlat);
    
    F_theta_coslat = F_theta.*coslat;
    dFdlat = ((Re*coslat).^-1).*[(F_theta_coslat(1,:,:)-F_theta_coslat(2,:,:))/dlat; (F_theta_coslat(1:MM-2,:,:)-F_theta_coslat(3:MM,:,:))/(2*dlat); (F_theta_coslat(MM-1,:,:)-F_theta_coslat(MM,:,:))/dlat];
    dFdp = [(F_p(:,2,:)-F_p(:,1,:))/(lev(2)-lev(1)), (F_p(:,3:LL,:)-F_p(:,1:LL-2,:))./repmat((lev(3:LL)-lev(1:LL-2))',[MM 1 TT]), (F_p(:,LL,:)-F_p(:,LL-1,:))/(lev(LL)-lev(LL-1))];
    divF = dFdlat + dFdp;
    
    if dd == DD
    	TT_use = TTT-((dd-1)*TT_old);

        Fp_surf((dd-1)*TT_old+1:TTT,1) = squeeze(sum(F_p(M_20:M_65,L_850,1:TT_use).*coslat(M_20:M_65,1,1:TT_use))/sum(coslat(M_20:M_65)));
        Fp_trop((dd-1)*TT_old+1:TTT,1) = squeeze(sum(F_p(M_20:M_65,L_200,1:TT_use).*coslat(M_20:M_65,1,1:TT_use))/sum(coslat(M_20:M_65)));
        divF_strat((dd-1)*TT_old+1:TTT,1) = squeeze(sum(divF(M_45:MM,L_10,1:TT_use).*coslat(M_45:MM,1,1:TT_use))/sum(coslat(M_45:MM)));

        Fp_surf((dd-1)*TT_old+1:TTT,2) = squeeze(sum(F_p(M_n65:M_n20,L_850,1:TT_use).*coslat(M_n65:M_n20,1,1:TT_use))/sum(coslat(M_n65:M_n20)));
        Fp_trop((dd-1)*TT_old+1:TTT,2) = squeeze(sum(F_p(M_n65:M_n20,L_200,1:TT_use).*coslat(M_n65:M_n20,1,1:TT_use))/sum(coslat(M_n65:M_n20)));
        divF_strat((dd-1)*TT_old+1:TTT,2) = squeeze(sum(divF(1:M_n45,L_10,1:TT_use).*coslat(1:M_n45,1,1:TT_use))/sum(coslat(1:M_n45)));
    else
        Fp_surf((dd-1)*TT+1:dd*TT,1) = squeeze(sum(F_p(M_20:M_65,L_850,:).*coslat(M_20:M_65,1,:))/sum(coslat(M_20:M_65)));
        Fp_trop((dd-1)*TT+1:dd*TT,1) = squeeze(sum(F_p(M_20:M_65,L_200,:).*coslat(M_20:M_65,1,:))/sum(coslat(M_20:M_65)));
        divF_strat((dd-1)*TT+1:dd*TT,1) = squeeze(sum(divF(M_45:MM,L_10,:).*coslat(M_45:MM,1,:))/sum(coslat(M_45:MM)));

    	Fp_surf((dd-1)*TT+1:dd*TT,2) = squeeze(sum(F_p(M_n65:M_n20,L_850,:).*coslat(M_n65:M_n20,1,:))/sum(coslat(M_n65:M_n20)));
    	Fp_trop((dd-1)*TT+1:dd*TT,2) = squeeze(sum(F_p(M_n65:M_n20,L_200,:).*coslat(M_n65:M_n20,1,:))/sum(coslat(M_n65:M_n20)));
        divF_strat((dd-1)*TT+1:dd*TT,2) = squeeze(sum(divF(1:M_n45,L_10,:).*coslat(1:M_n45,1,:))/sum(coslat(1:M_n45)));
    end
end

file_out = strcat(output_dir,'/EP_flux.nc');

nccreate(file_out,'Fp_surf','Dimensions',{'time',TTT,'hemisphere',2});
nccreate(file_out,'Fp_trop','Dimensions',{'time',TTT,'hemisphere',2});
nccreate(file_out,'divF_strat','Dimensions',{'time',TTT,'hemisphere',2});

ncwrite(file_out,'Fp_surf',Fp_surf);
ncwrite(file_out,'Fp_trop',Fp_trop);
ncwrite(file_out,'divF_strat',divF_strat);
