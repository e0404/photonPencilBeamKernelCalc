function ppbkc_visualizationScrtip(machineA,machineB)

if numel(machineA.data.kernel) ~= numel(machineB.data.kernel)
    error('different number of kernels in two files')
end

for i = 1:numel(machineA.data.kernel)

    if machineA.data.kernel(i).SSD ~= machineB.data.kernel(i).SSD
        error('SSDs do not match')
    end
    
    clf
    semilogx(machineA.data.kernelPos,machineA.data.kernel(i).kernel1,'b')
    hold on
    semilogx(machineB.data.kernelPos,machineB.data.kernel(i).kernel1,'b--')
    semilogx(machineA.data.kernelPos,machineA.data.kernel(i).kernel2,'r')
    semilogx(machineB.data.kernelPos,machineB.data.kernel(i).kernel2,'r--')
    semilogx(machineA.data.kernelPos,machineA.data.kernel(i).kernel3,'g')
    semilogx(machineB.data.kernelPos,machineB.data.kernel(i).kernel3,'g--')

    grid minor
    box on
    xlabel('[mm]')
    ylabel('a.u.')
    
    title(['SSD = ' num2str(machineA.data.kernel(i).SSD) 'mm'])
    
    legend({'kernel 1 A','kernel 1 B','kernel 2 A','kernel 2 B','kernel 3 A','kernel 3 B',})
    
    drawnow
    
    %pause(.01)
    
end
