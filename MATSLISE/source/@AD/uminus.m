function y = uminus(y)
% Unary minus for AD objects.
y.tc = -y.tc;
