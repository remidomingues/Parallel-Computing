function triangle(N,alpha)
Nf = floor(N/2);
% 
% 
% h = zeros(N,1);
% for f=1:Nf
%     P = 1/f^alpha * randn;
%     l = N/f;
%     for x=1:N
% %         h(x) = h(x) + P*(1 - 2*abs(mod(x/l+1/2,2) - 1));
% %         h(x) = (1 - 2*abs(mod(x/l+1/2,2) - 1));
%     end
%     scatter(1:N,h);
%     pause
% end

h = zeros(N,1);
h_old = zeros(N,1);
for f=1:Nf
    h_old = h;
    P = 1/f^alpha * randn;
    l = N/f/2;
%     hold on
    for x=1:N
%         h(x) = h(x) + P/l*(l - 2*abs(mod(x+1/2*l,2*l) - l));
        h(x) = (l - 2*abs(mod(x+1/2*l,2*l) - l));
    end
    scatter(1:N,h+h_old);
    hold on
    scatter(1:N,h-h_old);
    hold off
    pause
end

% 

% h = zeros(N,1);
% for f=1:Nf
%     l = (N/f);
%     h = zeros(N,1); %RM
%     P = (l/f)^alpha * randn;
%     P = 1; %RM
%     height = 0;
%     step = 2*P/l;
%     for x=0:ceil(l/2)
%         height = height+step;
%         sign = +1;
%         for i=0:f
%             left = round(i*l);
%             right = round(l*(i+1));
%             if (left~=0 || x~=0) && x+left < N
%                 h(left+x) =  h(left+x) + sign*height;
%             end
%             if (mod(l,2)~=0 || x~= ceil(l/2)) && (right-x > 0) && (right-x < N)
%                 h(right-x) = h(right-x) + sign*height;
%             end
%             sign = -sign;
%         end 
%     end
%     scatter(1:length(h),h);
% %     ylim([-1 1]);
%     xlim([0 N]);
%     pause
% end








end

