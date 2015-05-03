function lines(N, quality)
Nrays = floor(quality*N*N);
delta = 1;
h = zeros(N,N);


tic
for i=1:Nrays
    
    %generate random values
    y1 = rand; x1 = rand; x2 = rand; y2 = rand;
    m = (y2-y1)/(x2+x1);
    if rand>0.5, sign = 1;
    else sign = -1;
    end
    for x=1:N
        ylim = floor(m*(x - N*x1) + N*y1);
        if ylim>N, ylim = N;
        end
        for y=1:ylim
            h(x,y) = h(x,y) + sign*delta;
        end
    end
%     if mod(i,100)==1
%         imagesc(h)
%         colormap(copper);
%         pause
%     end
end
toc



imagesc(h)
colormap(copper);
% plotTerrain(h)

end
