function  plot_spectrum(signal, Fs)
%PLOT_SPECTRUM Summary of this function goes here
%   Detailed explanation goes here

    n = 2^nextpow2( length( signal ) );         % number of samples
    f = Fs * ( 0 : (n/2) ) / n;                 % frequency range                

    input_fft = fft( signal );
    input_fft_power = abs( input_fft / n );     % power of the DFT
    
    plot( f, input_fft_power( 1 : n/2 + 1 ) )
    xlabel('Frequency')
    ylabel('Power')

end

