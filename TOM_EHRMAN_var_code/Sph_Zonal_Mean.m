% This function generates a zonal mean remap matrix using spherical harmonic decomposition 
% given a latitude, output latitude, and order (L)

function [ZM, ZM_nat] = Sph_Zonal_Mean(lat,lat_out,L)

NN = length(lat); % # of input lats

LL = sum((1:L+1)); % sum of 1+2+3+...L+1

l = zeros(LL,1); % 0,1,1,2,2,2,3,3,3,3,...
m = zeros(LL,1); %0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,....
P = zeros(LL,NN);

sinlat = sind(lat);
for ll = 0:L
    l(sum((0:ll))+1:sum((0:ll+1))) = ll*ones(ll+1,1);
    m(sum((0:ll))+1:sum((0:ll+1))) = (0:ll);
    P(sum((0:ll))+1:sum((0:ll+1)),:) = legendre(ll,sinlat);
end

L0 = L+1;

Y0 = zeros(L0,NN);
l0 = 1;
for ll = 1:LL
    if m(ll) == 0
        coef = sqrt(((2*l(ll)+1)/(2*pi)));
        Y = coef*P(ll,:);
        Y0(l0,:) = Y(:);
        l0 = l0+1;
    end
end

NNN = length(lat_out);

P = zeros(LL,NNN);

sinlat = sind(lat_out);
for ll = 0:L
    P(sum((0:ll))+1:sum((0:ll+1)),:) = legendre(ll,sinlat);
end

Y1 = zeros(L0,NNN);
l0 = 1;
for ll = 1:LL
    if m(ll) == 0
        coef = sqrt(((2*l(ll)+1)/(2*pi)));
        Y = coef*P(ll,:);
        Y1(l0,:) = Y(:);
        l0 = l0+1;
    end
end

ZM = Y1'*(Y0'\eye(NN));
ZM_nat = Y0'*(Y0'\eye(NN));
