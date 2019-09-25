function plot_transfer_function_nofig( tf, ~, nunits )
%PLOT_TRANSFER_FUNCTION Plots bode of a transfer function with markers
%
%   tf                - The transfer function (created using tf)
%   frequency_markers - A matrix of frequencies in Hz
%
%   Example:
%       plot_transfer_function( tf([1000], [1 1000]), [10 1000 10000] );

    persistent Aw nui nui_max
    if ( isempty( Aw ) )

        nui_max = nunits + 1;
        Aw = cell(nui_max, 1);
        nui = 1;

    end

    x_space = logspace(1,5,5000); % 5000 points between 10^1 and 10^5
    x_space = 2 * pi * x_space; % to rad / sec
    [mag,~,wout] = bode(tf,x_space);
    mag = squeeze(mag);
    wout = squeeze(wout);
    mag = 20*log10(mag);
    wout = wout/2/pi;
    semilogx(wout,mag,'-');
    axis([min(wout) max(wout) min(mag)-30 max(mag)+30]);

    xlabel('Frequency (Hz)', 'FontSize', 18);
    ylabel('Magnitude (dB)', 'FontSize', 18);
    grid on;
    hold on;

    if ( nui < nui_max )
        Aw{nui} = ['|T' num2str(nui) '(jω)|'];
        nui = nui + 1;
    else
        Aw{nui} = 'Total TF: |T(jω)|';
        legend(Aw,'Location','best','FontSize',12);
        set(gca,'FontSize',14);
    end

end
