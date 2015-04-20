function FourierFractal(N,alpha)
Nf = ceil(N/4);
% Px = zeros(1,Nf);
% Py = zeros(1,Nf);
h = zeros(N,N);
P = zeros(Nf,Nf);

for fx=1:Nf
    for fy=1:Nf
        P(fx,fy) = N^2 * (fx^2+fy^2)^(-alpha/2) * randn * exp(2i * pi * rand); 
%         P(fx,fy) = (fx^2+fy^2)^(-alpha/2) * (randn + 1i*randn);
    end
end
P(1,1) = 0;
% for f=2:Nhalf
%     Px(f) = f^(-alpha) * randn * exp(2i * pi * rand);
%     Py(f) = f^(-alpha) * randn * exp(2i * pi * rand);
% end
% 
% P = Px'*Py;
% h = ifft(Px)' * ifft(Py);

h = ifft2(P);

% for x=1:N
%     for y=1:N
%         for fx=1:Nhalf
%             for fy=1:Nhalf
% %                 h(x,y) = h(x,y) + 1/(sqrt(fx^2+fy^2)^alpha) * real(P(fx,fy) * exp(2i*pi*(fx*x + fy*y)/N));
%                 h(x,y) = h(x,y) + 1/(sqrt(fx^2+fy^2)^alpha) * real(Px(fx) * Py(fy) * exp(2i*pi*(fx*x + fy*y)/N));
%             end
%         end
%     end
% end
% figure
% abs(P)



plotTerrain(real(h));
end
    
