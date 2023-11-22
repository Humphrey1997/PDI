fsize = 15;

subplot(4,4,1)
plot(S_ns,real(S_eps))
title('real-eps')
set(gca,'fontsize',fsize)

subplot(4,4,2)
plot(S_ns,imag(S_eps))
title('imag-eps')
set(gca,'fontsize',fsize)

subplot(4,4,3)
plot(S_ns,real(S_ep1))
title('real-ep1')
set(gca,'fontsize',fsize)

subplot(4,4,4)
plot(S_ns,imag(S_ep1))
title('imag-ep1')
set(gca,'fontsize',fsize)

subplot(4,4,5)
plot(S_ns,real(S_QLs))
title('real-QLs')
set(gca,'fontsize',fsize)

subplot(4,4,6)
plot(S_ns,imag(S_QLs))
title('imag-QLs')
set(gca,'fontsize',fsize)

subplot(4,4,7)
plot(S_ns,real(S_QL1))
title('real-QL1')
set(gca,'fontsize',fsize)

subplot(4,4,8)
plot(S_ns,imag(S_QL1))
title('imag-QL1')
set(gca,'fontsize',fsize)

subplot(4,4,9)
plot(S_ns,real(S_Lse))
title('real-Lse')
set(gca,'fontsize',fsize)

subplot(4,4,10)
plot(S_ns,imag(S_Lse))
title('imag-Lse')
set(gca,'fontsize',fsize)

subplot(4,4,11)
plot(S_ns,real(S_Lsi))
title('real-Lsi')
set(gca,'fontsize',fsize)

subplot(4,4,12)
plot(S_ns,imag(S_Lsi))
title('imag-Lsi')
set(gca,'fontsize',fsize)

subplot(4,4,13)
plot(S_ns,real(S_miu))
title('real-miu')
set(gca,'fontsize',fsize)

subplot(4,4,14)
plot(S_ns,imag(S_miu))
title('imag-miu')
set(gca,'fontsize',fsize)

subplot(4,4,15)
plot(S_ns,real(S_ep))
title('real-ep')
set(gca,'fontsize',fsize)

subplot(4,4,16)
plot(S_ns,imag(S_ep))
title('imag-ep')
set(gca,'fontsize',fsize)