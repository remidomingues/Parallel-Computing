function FourierFractal(N,alpha)
Nhalf = ceil(N/2);
Px = zeros(1,Nhalf);
Py = zeros(1,Nhalf);
h = zeros(N,N);
P = zeros(Nhalf,Nhalf);

P = randn(Nhalf,Nhalf) +i*randn(Nhalf,Nhalf);

for f=1:Nhalf
    Px(f) = f^(-alpha) * randn * exp(2i * pi * rand);
    Py(f) = f^(-alpha) * randn * exp(2i * pi * rand);
end
% 
 P = Px'*Py;

h = ifft(Px)' * ifft(Py);

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
figure
surf(real(h))
% figure
% surf(real(f))


% plotTerrain(h);
end
    
