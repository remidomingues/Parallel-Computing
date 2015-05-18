function demo(N)
midpoint = zeros(N,1);
lines = zeros(N,1);
fourier = zeros(N,1);

k=0;
damping = 0.9;

figure
for iter=0:N/2
    % midpoint
    k = 2^iter;
    if k < (N-1)
        d = (N-1)/(k);
        for i=(d/2+1):d:(N-1)
            midpoint(i) = (midpoint(i-d/2) + midpoint(i+d/2))/2 + ...
                1/k^damping*(2*rand-1);
        end
    end
    
    
    % line
    lim = floor(N*rand);
    if rand>0.5, sign = 1;
    else sign = -1;
    end
    if rand<0.5
        lines(1:lim) = lines(1:lim) + sign;
    else
        lines(lim:N) = lines(lim:N) + sign;
    
    % fourier
    f = iter + 1;
    p =  1/f  * exp(2i*pi*rand);
    for x=1:N
        fourier(x) = fourier(x) + real(p * exp(2i*pi*f*x/N));
    end

%plots
subplot(3,1,1)
scatter(1:N, midpoint,'b.');
title('Midpoint');
xlim([1 N])
ylim([-1 1]);

subplot(3,1,2)
scatter(1:N, lines, 'k.');
title('Linear displacement');
xlim([1 N])
m = max(abs(lines));
ylim([-m m]);

subplot(3,1,3)
plot(1:N, fourier);
title('Fourier transform');
xlim([1 N])
ylim([-2 2]);

pause

end
end


	
