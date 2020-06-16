
settings = struct();
for i_part = 1:length(parts)
    for i_chan = 1:n_chans
            fooof_eeg(i_sub,i_chan,i_cond) = fooof(F, bosc_spectra_baseline(:,i_sub,i_chan,i_cond), [0.1,30], settings,'1');
            fooof_eeg_spectra(:,i_sub,i_chan,i_cond) = fooof_eeg(i_sub,i_chan,i_cond).power_spectrum;
            fooof_eeg_fooofed_spectra(:,i_sub,i_chan,i_cond) = fooof_eeg(i_sub,i_chan,i_cond).fooofed_spectrum;
            fooof_eeg_bg_spectra(:,i_sub,i_chan,i_cond) = fooof_eeg(i_sub,i_chan,i_cond).bg_fit;
    end
end

