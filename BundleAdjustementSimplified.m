function [Rwi_c, Twi_c, Uw] = BundleAdjustementSimplified(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, maxIt, nCam)

lambda = 1e-3;
lambdaMin = 1e-5;
lambdaMax = 1e5;
%compute cost
c = compute_cost(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam);

fprintf('Iter:\t Error:\t lambda:\n');
fprintf('%d\t%f\t%f\n',0,c,lambda);

for i = 1:maxIt
    
    %compute Jacobian and residuals
    r = compute_residuals(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam);
    J = compute_Jacobian(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam);
    
    %update parameters
    [Rwi_c, Twi_c, Uw, c, lambda, success] = updateParameters(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, r, J, lambda, lambdaMax, lambdaMin, c, i, nCam);
    
    if(~success)
        break;
    end

end

end

function c = compute_cost(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam)

r = compute_residuals(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam);
c = r'*r;

end

function r = compute_residuals(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam)


r = [];

p_vis = p{nCam}(:,tracks_cur{1}.point2D_ids);
    
    Ui = Rwi_c'*(Uw(:,tracks_cur{1}.point3D_ids) - Twi_c);
    p_vis_pred_hom = K*Ui./repmat(Ui(3,:),3,1);
    p_vis_pred = p_vis_pred_hom(1:2,:);
    
    r = [r; p_vis(:) - p_vis_pred(:)];

    
end


function J = compute_Jacobian(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, nCam)

nb_Pt_2D = length(tracks_cur{1}.point3D_ids);

J = zeros(2*nb_Pt_2D, 6*nCam+3*nTracks);
    
    
[G1, G2, G3] = generators_SO3();
 
l=1;

Ui = Rwi_c'*(Uw(:,tracks_cur{1}.point3D_ids) - Twi_c);
    
    

for t = 1 :nb_Pt_2D %for each 2D point

    A = [1/Ui(3,t) 0 -Ui(1,t)/(Ui(3,t)^2); 0 1/Ui(3,t) -Ui(2,t)/(Ui(3,t)^2); 0 0 0];

    %3D point derivative
    J(l:l+1, 1+ 3*(tracks_cur{1}.point3D_ids(t)-1):3*(tracks_cur{1}.point3D_ids(t))) = K(1:2,:)*A*Rwi_c';

    %rotation derivative
    J(l:l+1,3*nTracks+1:3*nTracks+3) = [K(1:2,:)*A*G1'*Ui(:,t), K(1:2,:)*A*G2'*Ui(:,t), K(1:2,:)*A*G3'*Ui(:,t)];

    %translation derivative
    J(l:l+1,3*nTracks+4:3*nTracks+6) = -K(1:2,:)*A*Rwi_c';

    l = l+2;
end


end

function [G1, G2, G3] = generators_SO3()

G1 = [0  0 0; 0 0 -1;  0 1 0];
G2 = [0  0 1; 0 0  0; -1 0 0];
G3 = [0 -1 0; 1 0  0;  0 0 0];

end

function w_hat = HatSO3(w)

w_hat = [0 -w(3) w(2);...
    w(3) 0 -w(1);...
    -w(2) w(1) 0];

end

function [Rwi_new_c, Twi_new_c, Uw_new, c_new, lambda, success] = updateParameters(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, r, J, lambda, lambdaMax, lambdaMin, c_prev, iter, nCam)

success = false;


while(lambda<lambdaMax)
    
    %solve linear system
    delta = mldivide(J'*J + lambda*eye(3*nTracks + 6*nCam), J'*r);
    
    
    %update variables
    % Uw_new = Uw + reshape(delta(1:3*nTracks),3,nTracks);
    Uw_new = Uw; 
    
    Rwi_new_c = Rwi_c*expm(HatSO3(delta(3*nTracks+1:3*nTracks+3)));
    Twi_new_c = Twi_c + delta(3*nTracks+4:3*nTracks+6);
    
    
    %compute cost
    c_new = compute_cost(Rwi_new_c, Twi_new_c, Uw_new, tracks_cur, K, p, nTracks, nCam);
    fprintf('%d\t%f\t%f\n',iter,c_new,lambda);

    if((c_new+1e-5) < c_prev)
        success = true;
        if(lambda>lambdaMin)
            lambda = lambda/2;
        end
        break;
    else
        lambda = lambda*2;
    end
end

if(~success)
    Rwi_new_c = Rwi_c;
    Twi_new_c = Twi_c;
    Uw_new = Uw;
    c_new = c_prev;
end

end




