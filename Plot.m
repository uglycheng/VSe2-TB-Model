for Vdda=-2:2
    for Vddb=-2:2
        for Vddc=-2:2
            for Vpda=-2:2
                for Vpdb=-2:2
                    name=strcat('V',num2str(Vdda),'_',num2str(Vddb),'_',num2str(Vddc),'_',num2str(Vpda),'_',num2str(Vpdb));
                    run('TBH');
                    scatter(x_cordi,Band,7,'filled');
                    print('-djpeg','-r100',name);
                end
            end
        end
    end
end