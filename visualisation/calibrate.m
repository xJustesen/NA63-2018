function [minimum,val, k, f, dat, sim, dat_err] = calibrate(E, amorphsim, bgsim, amorphdat, bgdat, NsimAmorph, NsimBg, Ndat1, Ndat2)
    amorph_sim = hist(amorphsim(amorphsim < max(E) & amorphsim > min(E)), E);
    bg_sim = hist(bgsim(bgsim < max(E) & bgsim > min(E)), E);
    sim = E .* (amorph_sim/NsimAmorph - bg_sim/NsimBg);

    amorph_dat = hist(amorphdat(amorphdat < max(E) & amorphdat > min(E)), E);
    bg_dat = hist(bgdat(bgdat < max(E) & bgdat > min(E)), E);
    dat = E .* (amorph_dat/Ndat1 - bg_dat/Ndat2);
    dat_err = E .* sqrt(amorph_dat/Ndat1^2 + bg_dat/Ndat2^2);
   
    if (E(end) == 20)
        limit = -1.5e-4;
    elseif (E(end) == 40)
        limit = -4e-4;
    else
        limit = -6e-4;
    end
%     limit = -1e9;
    chi2 = @(eff) norm(dat(dat > limit) - eff .* sim(dat > limit));
%         chi2 = @(eff) trapz(dat(dat > limit) - eff .* sim(dat > limit));
    
    f = zeros(100,1);
    k = linspace(0,2,100);
    for i = 1:100
        f(i) = chi2(k(i));
    end

    minimum = fminsearch(chi2, 0.5); % fast, local
%     minimum = fzero(chi2, 0.5);
    val = chi2(minimum);
    
%     figure
%     hold on
%     plot(E, E .* amorph_sim/NsimAmorph,'b-','linewidth',1.5)
%     plot(E, E .* bg_sim/NsimBg,'r--','linewidth',1.5)
%     plot(E, sim,'g-.','linewidth',1.5)
%     plot(E, E .* amorph_dat/Ndat1,'bo','MarkerFaceColor','b')
%     plot(E, E .* bg_dat/Ndat2,'rs','MarkerFaceColor','r')
%     errorbar(E, dat,dat_err,'kh','MarkerFaceColor','g')
%     plot(E, limit*ones(size(E)),'k-','linewidth',1.5)
%     legend('amorf+bg sim','bg sim','amorph-bg sim','amorph+bg dat','bg dat','amorph - bg dat')
%     box on; grid on;
    
end
