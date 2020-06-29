%% SHORT SCRIPT TO NUMERICALLY SOLVE ICE FLOW/SLIDING EQUATIONS AND INVERT FOR ICE THICKNESS


%% inputs

n = 3;  %Glen's flow constant

A_c = 2.5e-24;  %Arrhenius creep constant

f = 0.9; %lateral drag correction ('shape factor')

p_i = 917; %density of ice

g = 9.79; %gravity

spery = 60*60*24*365; %seconds per year

C1_1 = 0.000/spery; %0.0012 Sliding constant

C2_1 = 0; %Basal water fraction

% threshold_1 = 6000; %if threshold added
% 
% C1_2 = 0.0005/spery; %0.0012
% 
% C2_2 = 0.1;
% 
% threshold_2 = 5500;
% 
% C1_3 = 0.0010/spery; %0.0012
% 
% C2_3 = 0.25;


% %Smooth slope along 150m coupling length
% slope = movmean(slope,3,2);
% slope = movmean(slope,3,1);



% %Make slope clip same size as velocity matrix if necessary
% Slope_clip = interp2(Slope_clip, linspace(1, size(Slope_clip,2), size(Ice_velocity_resamples,2)).', linspace(1, size(Slope_clip,1), size(Ice_velocity_resamples,1)));


% Reduce size to speed up calculations if necessary
% slope = interp2(slope, linspace(1, size(slope,2), size(slope,2)/5).', linspace(1, size(slope,1), size(slope,1)/5));
% velocity = interp2(velocity, linspace(1, size(velocity,2), size(velocity,2)/5).', linspace(1, size(velocity,1), size(velocity,1)/5));

% save_volume = NaN(5,5); %If iterating
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
        
%         if DEM_long(loop,1) > threshold_1  %If threshold added
            syms H

        S = vpasolve((H^(n+1)*2*A_c*(f*p_i*g*sin(slope_long(loop,1)))^n/(n+1))...
     + (H^2*sin(slope_long(loop,1))^2*((C1_1*(f*p_i*g)^2)/((1-C2_1)*(H*p_i*g))))...
     == velocity_long(loop,1), H);
 
         storage_mat(loop,1) = S(2,1);
         
%         elseif DEM_long(loop,1) <= threshold_1 && DEM_long(loop,1) > threshold_2
%         syms H
% 
%          S = vpasolve((H^(n+1)*2*A_c*(f*p_i*g*sin(slope_long(loop,1)))^n/(n+1))...
%              + (H^2*sin(slope_long(loop,1))^2*((C1_2*(f*p_i*g)^2)/((1-C2_2)*(H*p_i*g))))...
%        == velocity_long(loop,1), H);
%  
%  storage_mat(loop,1) = S(2,1);
%  
%         elseif DEM_long(loop,1) <= threshold_2
%         syms H
% 
%          S = vpasolve((H^(n+1)*2*A_c*(f*p_i*g*sin(slope_long(loop,1)))^n/(n+1))...
%              + (H^2*sin(slope_long(loop,1))^2*((C1_3*(f*p_i*g)^2)/((1-C2_3)*(H*p_i*g))))...
%        == velocity_long(loop,1), H);
%  
%  storage_mat(loop,1) = S(2,1);
%         end
%              
    end
end
 toc

inverted_thickness=reshape(storage_mat,size(slope));




