fsize = 15;

subplot(2,2,1)
plot(S_ns,real(S_eps))
title('real-eps')
set(gca,'fontsize',fsize)

subplot(2,2,2)
plot(S_ns,imag(S_eps))
title('imag-eps')
set(gca,'fontsize',fsize)

subplot(2,2,3)
plot(S_ns,real(S_ep1))
title('real-ep1')
set(gca,'fontsize',fsize)

subplot(2,2,4)
plot(S_ns,imag(S_ep1))
title('imag-ep1')
set(gca,'fontsize',fsize)


