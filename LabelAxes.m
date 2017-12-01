function handl = LabelAxes(h,numb,FontSize,xp,yp)
% This function puts a letter label in the top left corner of the
% specified axes

if (exist('h')==0)
    h = gca;
end
if (exist('numb')==0)
    numb = 1;
end
if (exist('xp')==0)
    xp = 0.02;
end
if (exist('yp')==0)
    xp = 0.94;
end


txt = {'(a)','(b)','(c)','(d)','(e)','(f)',...
       '(g)','(h)','(i)','(j)','(k)','(l)',...
       '(m)','(n)','(o)','(p)','(q)','(r)',...
       '(s)','(t)','(u)','(v)','(w)','(x)'};
xlims = get(h,'xlim');xL = xlims(2)-xlims(1);
ylims = get(h,'ylim');yL = ylims(2)-ylims(1);

if (exist('FontSize')==0)
    text(xlims(1)+xL*xp,ylims(1)+yL*yp,txt(numb),'BackgroundColor','w','margin',0.5);
else    
    text(xlims(1)+xL*xp,ylims(1)+yL*yp,txt(numb),'BackgroundColor','w','margin',0.5,'FontSize',FontSize);
end

end
