function [Rwi_BA_c,Twi_BA_c] = getPosiCam(Rwi_BA_c,Twi_BA_c,Uw_BA, tracks_cur, tracks, nCam,K,p)
    
    % On recherche les matches 
    [Li1, Loc2] = ismember(tracks_cur{nCam-1}.oldpoint3D_ids, tracks{nCam}.point3D_ids);
    
    nTracks=sum(Li1);
    
    tracks_cur_nCam = tracks{nCam};
    tracks_cur_nCam.point2D_ids = tracks{nCam}.point2D_ids(Loc2(Li1));
    tracks_cur_nCam.point3D_ids=1:nTracks;
       
    % On cherche les correspondnace entre les points 3D qio sont dans  im2 et im3 afin
    % de faire la pojection et estimer R et T 
    
    Uw_commun= Uw_BA(:,tracks_cur{nCam - 1}.point3D_ids(Li1));
    
    maxIt = 20;
    [RwiCam, TwiCam] = BundleAdjustementSimplified(Rwi_BA_c{nCam-1}, Twi_BA_c{nCam-1}, Uw_commun, {tracks_cur_nCam}, K, p, nTracks, maxIt, nCam);
    Rwi_BA_c = [Rwi_BA_c RwiCam]; 
    Twi_BA_c = [Twi_BA_c TwiCam];
end



