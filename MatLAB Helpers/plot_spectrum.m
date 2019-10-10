function  plot_spectrum(signal, Fs)
%PLOT_SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

    n = length( signal );    
    f = 0 : Fs / n : Fs - Fs / n;       % frequency range

    signal_fft = fft( signal );    
    singal_fft_power = 2 * abs( signal_fft / n);     % power of the DFT
    
    plot( f, singal_fft_power )
    xlim( [0 Fs/2] );
    xlabel('Frequency (Hz)')
    ylabel('|FT(f)| (V)')

end

