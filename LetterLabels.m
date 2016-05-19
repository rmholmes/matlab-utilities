function handl = LetterLabels(n,FontSize)
% This function puts letter label text in a figure for labels
% different subplots.

if (exist('FontSize')==0)
    FontSize = 20;
end

txt = {'(a)','(b)','(c)','(d)','(e)','(f)',...
       '(g)','(h)','(i)','(j)','(k)','(l)'};
handl = zeros(n,1);
for ii = 1:n
    %handl(ii) = uicontrol('style','text');
    set(handl(ii),'String',txt{ii},'FontSize',FontSize,'BackgroundColor', ...
                'w','Position',[10 10 35 30]);
end
end
