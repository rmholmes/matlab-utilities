function PDF = FloatPDF(Fval,BinLim,NoLvls,FilterWindow);
%PDF: This function calculates the probability density function of
%a set of floats with respect to float number. Fval is a vector
%containing the values of the PDF variable on each float, and
%BinLim = [mnval mxval] specifies the maximum and minimum values
%for the x-axis.

%Parameters for pdf:
% $$$ NoLvls = 750;
% $$$ FilterWindow = 5; %must be odd!
OneSide = (FilterWindow-1)/2;
mnval = BinLim(1);
mxval = BinLim(2);

%Construct cdf:
dval = (mxval-mnval)/NoLvls;
CDF = zeros(NoLvls+1,2);
CDF(:,1) = mnval:dval:mxval;
for i=2:(NoLvls+1)
    CDF(i,2) = CDF(i-1,2)+length(find(Fval(Fval<=CDF(i,1) & Fval>CDF(i-1,1))));
end

%take derivative for pdf:
PDFraw = zeros(NoLvls,2);
PDFraw(:,1) = (CDF(2:end,1)+CDF(1:(end-1),1))/2;
PDFraw(:,2) = (CDF(2:end,2)-CDF(1:(end-1),2))./((CDF(2:end,1)- ...
                                                 CDF(1:(end-1),1)));
%Filter:
PDF = zeros(NoLvls-2*OneSide,2);
PDF(:,1) = PDFraw((OneSide+1):(NoLvls-OneSide),1);
for i = (OneSide+1):(NoLvls-OneSide)
    PDF(i-OneSide,2) = sum(PDFraw((i-OneSide):(i+OneSide),2))/ FilterWindow;
end

%Normalize:
Total = sum((PDF(1:(end-1),2)+PDF(2:end,2))/2.*(PDF(2:end,1)-PDF(1: ...
                                                  (end-1),1)));
PDF(:,2) = PDF(:,2)/Total;

end
