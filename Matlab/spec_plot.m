function [outputArg1,outputArg2] = spec_plot(Spec,N,edge,demod,o,DC)
%SPEC_PLOT Summary of this function goes here
%   Detailed explanation goes here
    
    figure
    Freq = 1:N;
    imagesc([1:length(Spec)],Freq,abs(Spec));
    if(edge)
        JET = colormap(gray);
        hold on
        if(demod)
            for temp = [N/8:N/8:N-(N/8)]
                plot(temp.*ones(1,N),[1:N],'linewidth',2);
            end
        else
            if(DC)
                plot(o(:),ones(1,2*size(o,1)),'o','MarkerSize',30,'linewidth',2);
                plot(N/2.*ones(1,length(Spec)),'linewidth',5);
            else
                plot(o(:)',ones(1,size(o,1)),'o','MarkerSize',30,'linewidth',2);
            end
        end
    else
        hold on
        JET = colormap(jet);
        if(~demod)
            plot(o(:)',ones(1,size(o,1)),'o','MarkerSize',30,'linewidth',2);
        end
    end
    
    mycolor = cat(1, [1 1 1], JET(2:end,:));
%     colormap (copper);   % Set the current colormap Style 
    colormap (JET);
    box on;
    colorbar off;
    set(gcf,'position',[846.6,340.2,414.4,364]);
    set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
%     ylim([Freq(1),Freq(end)]);
    ylim([Freq(1),Freq(end)]);
    view(0,90);
    set(gca,'YDir','normal');
    title('Chirplet TF of Current window','FontSize',30);
    xlabel('Sample','FontSize',30);
%     ylabel('Freq. / Hz','FontSize',30);
    ylabel('Tx no.','FontSize',30);
end

