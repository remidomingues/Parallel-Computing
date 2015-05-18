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
%     if (i==1 || mod(log10(i),1)==0)
%         imagesc(h)
%         colormap(copper);
%         titlestr = sprintf('Iterations: %d',i);
%         title (titlestr);
%         axis off;
%         axis equal;
%         pause
%     end
end
toc

%  plotTerrain(h);

end
