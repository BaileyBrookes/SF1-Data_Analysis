subplot(4,2,1)
plot(t,S)
xlabel('Time (s)')
ylabel('Amplitude')
title('Normal signal-Time Domain')
subplot(4,2,2)
plot(f_norm, data{1})
ylabel('Normalised Amplitude')
xlabel('Normalsied frequency (rads^{-1})')
title('Normal signal-Frequency Domain')

subplot(4,2,3)
plot(t,X1)
xlabel('Time (s)')
ylabel('Amplitude')
title('Random Modulation -Time Domain')
subplot(4,2,4)
plot(f_norm, data{2})
ylabel('Normalised Amplitude')
xlabel('Normalsied frequency (rads^{-1})')
title('Random Modulation-Frequency Domain')

subplot(4,2,5)
plot(t,X2)
xlabel('Time (s)')
ylabel('Amplitude')
title('Linear Increase-Time Domain')
subplot(4,2,6)
plot(f_norm, data{3})
ylabel('Normalised Amplitude')
xlabel('Normalsied frequency (rads^{-1})')
title('Linear Increase-Frequency Domain')

subplot(4,2,7)
plot(t,X3)
xlabel('Time (s)')
ylabel('Amplitude')
title('Periodic Modulation-Time Domain')
subplot(4,2,8)
plot(f_norm, data{4})
ylabel('Normalised Amplitude')
xlabel('Normalsied frequency (rads^{-1})')
title('Periodic Modulation-Frequency Domain')

set(findall(gcf,'-property','FontSize'),'FontSize',8)
