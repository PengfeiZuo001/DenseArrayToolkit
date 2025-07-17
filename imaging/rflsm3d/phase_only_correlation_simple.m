function [delay_time, poc_result, lag] = phase_only_correlation_simple(sig1, sig2, dt)
    % 保证输入为列向量
    sig1 = sig1(:); sig2 = sig2(:);
    N = min(length(sig1), length(sig2));
    sig1 = sig1(1:N); sig2 = sig2(1:N);
    sig1 = sig1 - mean(sig1); sig2 = sig2 - mean(sig2);
    win = hanning(N); % 列向量
    sig1 = sig1 .* win; sig2 = sig2 .* win;
    F1 = fft(sig1); F2 = fft(sig2);
    cross_power = F2 .* conj(F1); % 保持符号一致
    R = cross_power ./ (abs(cross_power) + 1e-12);
    poc = real(ifft(R));
    poc = fftshift(poc); % 列向量
    lag = (-floor(N/2)):(N - floor(N/2) - 1); lag = lag(:); % 列向量
    [~, max_idx] = max(poc);
    delay_time = lag(max_idx) * dt;
    poc_result = poc;
end