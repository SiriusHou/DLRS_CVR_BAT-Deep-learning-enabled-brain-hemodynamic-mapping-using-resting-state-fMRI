function [minRow,minCol] = zach_getMinCoords(resid)

[M,I] = min(resid);

[M2,I2] = min(M);

minRow = I(I2);

minCol = I2;