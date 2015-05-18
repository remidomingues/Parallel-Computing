function FourierFractal(N,alpha)
Nf = ceil(N/4);
h = zeros(N,N);
P = zeros(Nf,Nf);

for fx=1:Nf
    for fy=1:Nf
        P(fx,fy) = (fx^2+fy^2)^(-alpha) * randn * exp(2i * pi * rand);
    end
end
P(1,1) = 0;
h = fft2(P);

plotTerrain(real(h));
end
    
