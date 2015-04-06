function cosineFractal(N,alpha)
% f = randn(N,N);
h = zeros(N,N);

fx = randn(1,N);
fy = randn(1,N);

for x=1:N
    for y=1:N
        f(x,y) = (x^2+y^2)^(-alpha/2) * randn;
    end
end
figure
surf(log(abs(f)));
% f= fx' * fy;

h = dct2(f);
% for x=1:N
%     for y=1:N
%         for fx=1:N
%             for fy=1:N
%                 h(x,y) = h(x,y) + ...
%                     (fx^2+fy^2)^(-alpha/2) * f(fx,fy) * cos(pi*(fx+0.5)*x/N) * cos(pi*(fy+0.5)*y/N);
%             end
%         end
%     end
% end
figure
plotTerrain(h);
