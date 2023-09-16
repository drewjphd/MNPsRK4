tic; [M,t]=BrownSRK4; toc;
figure; plot(t,M); title('SRK4'); set(gca,'FontSize',16);

tic; [M,t]=BrownV2v2; toc;
figure; plot(t,M); title('Euler-Marayuma (Dan)'); set(gca,'FontSize',16); 