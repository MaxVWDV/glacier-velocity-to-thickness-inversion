%% INVERSION WITH REGELATION


%% inputs

n = 3;

A_c = 2.5e-24;

f = 0.9;

p_i = 917;

g = 9.79;

k_r = 2;

B_r = 7.4e-8;

L_f = 3.34e5;

R_r = 0.01; 


spery = 60*60*24*365;

% %Smooth slope along 150m coupling length
% slope = movmean(slope,3,2);
% slope = movmean(slope,3,1);



% %Make slope clip same size as velocity matrix if needed
% Slope_clip = interp2(Slope_clip, linspace(1, size(Slope_clip,2), size(Ice_velocity_resamples,2)).', linspace(1, size(Slope_clip,1), size(Ice_velocity_resamples,1)));


% Reduce size to speed up calculations if needed
% slope = interp2(slope, linspace(1, size(slope,2), size(slope,2)/5).', linspace(1, size(slope,1), size(slope,1)/5));
% velocity = interp2(velocity, linspace(1, size(velocity,2), size(velocity,2)/5).', linspace(1, size(velocity,1), size(velocity,1)/5));

% save_volume = NaN(5,5); if iterating
% save_summit = NaN(5,5);



% % iterate C1 and C2
% C1list = [0.0001,0.0005,0.001,0.002];
% for C1loop = 1:length(C1list)
% C2list = [0,0.1,0.25,0.5,0.95];
% for C2loop = 1:length(C2list)
    
%% Solve

slope_long=deg2rad(reshape(slope,[size(slope,1)*size(slope,2),1]));
velocity_long=reshape(velocity,[size(velocity,1)*size(velocity,2),1])/spery;
DEM_long=reshape(DEM,[size(DEM,1)*size(DEM,2),1]);

storage_mat = NaN(size(velocity_long,1),1);
% curr_pix = zeros(size(velocity_long,1)*4,1);


tic

for loop = 1:size(velocity_long,1)
    if ~isnan(velocity_long(loop,1))

            syms H

        S = vpasolve((H^(n+1)*2*A_c*(f*p_i*g*sin(slope_long(loop,1)))^n/(n+1))...
     + (H^((n+1)/2)*sin(slope_long(loop,1))^((n+1)/2)*((f*p_i*g)^((n+1)/2)*(2^((3-n)/2)/3^((n+1)/4))*((k_r*B_r*A_c)/(p_i*L_f))^0.5*(1/R_r)^(n+1)))...
     == velocity_long(loop,1), H);
 
         storage_mat(loop,1) = S(2,1);
         

             
    end
end
 toc

inverted_thickness=reshape(storage_mat,size(slope));

nanmean(inverted_thickness,'all')


figure;surf(inverted_thickness)






