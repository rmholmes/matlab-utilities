
%This script draws a line between two points on a matlab .fig file
%where the two points are p1 and p2 taken from a datatip output.

%slope:
m_l2p = (p2.Position(2)-p1.Position(2))/(p2.Position(1)-p1.Position(1));

%intercept:
b_l2p = p1.Position(2)-m_l2p*p1.Position(1);

%plot line:
hold on;
xlims_l2p = get(gca,'xlim');
xvec_l2p = xlims_l2p(1):((xlims_l2p(2)-xlims_l2p(1))/500):xlims_l2p(2);
plot(xvec_l2p,m_l2p*xvec_l2p+b_l2p,'-k');

'Intercept'
b_l2p
'Slope'
m_l2p