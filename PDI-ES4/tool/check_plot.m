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
plot(S_ns,real(S_NLs))
title('real-NLs')
set(gca,'fontsize',fsize)

subplot(4,4,10)
plot(S_ns,imag(S_NLs))
title('imag-NLs')
set(gca,'fontsize',fsize)

subplot(4,4,11)
plot(S_ns,real(S_NL1))
title('real-NL1')
set(gca,'fontsize',fsize)

subplot(4,4,12)
plot(S_ns,imag(S_NL1))
title('imag-NL1')
set(gca,'fontsize',fsize)

subplot(4,4,13)
plot(S_ns,real(S_NNLs))
title('real-NNLs')
set(gca,'fontsize',fsize)

subplot(4,4,14)
plot(S_ns,imag(S_NNLs))
title('imag-NNLs')
set(gca,'fontsize',fsize)

subplot(4,4,15)
plot(S_ns,real(S_NNL1))
title('real-NNL1')
set(gca,'fontsize',fsize)

subplot(4,4,16)
plot(S_ns,imag(S_NNL1))
title('imag-NNL1')
set(gca,'fontsize',fsize)

