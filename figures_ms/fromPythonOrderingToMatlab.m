% A function to use the Python Ordering, then switch to be consistent with Matlab's ordering
function newOrder = fromPythonOrderingToMatlab(originalInds,pythonOrdering,volume)
% This scripts has to do a bit of fiddling around -- you have to recreate the ordering in a volume
% to find the transformation between python and matlab's reshape, then apply it

temp = 1:length(originalInds);
gm2 = 0*volume;
gm2(originalInds) = temp;
% Annoying little change: Due to the way that matlab does it
mm = permute(gm2,[2 1 3]);
newInds = mm(find(mm));
newOrder = newInds(pythonOrdering);

