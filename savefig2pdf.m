function savefig2pdf(filename)
if nargin<1
    [FileName,PathName] = uiputfile;
    if FileName~=0
        filename = [PathName FileName];
    end
end
fig2save = gcf;
set(fig2save,'Units','Inches');
pos = get(fig2save,'Position');
set(fig2save,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if ~contains(filename,'.pdf')
    savefig([filename '.fig']);    
else
    print(fig2save,filename,'-dpdf','-r0')
end
end
