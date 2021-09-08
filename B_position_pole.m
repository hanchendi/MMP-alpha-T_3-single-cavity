figure; 
subplot(211)
plot(V2R, P0, 'r');
%xlim([.48 .51])
subplot(212)
plot(V2R, ee, 'ko');
%xlim([.48 .51])


figure;
plot(X, Y, 'k--.')
hold on
plot(real(zRm), imag(zRm), 'ro')
plot(real(zRl), imag(zRl), 'b*')
plot(real(z0), imag(z0), 'mp','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10)
%plot([real(Z); real(Zf)], [imag(Z); imag(Zf)], 'm')
axis equal; 