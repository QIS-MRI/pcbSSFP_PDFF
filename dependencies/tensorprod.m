function tprod = tensorprod(A,B, dimA, dimB)

szA = size(A); szDimA = szA(dimA);
dimsA = 1:1:length(szA);
szA(dimA) = [];
dimsA(dimA) = [];

permuteA = permute(A, [dimsA, dimA]);
permuted_reshapedA = reshape(permuteA, prod(szA), szDimA);

szB = size(B); szDimB = szB(dimB);
dimsB = 1:1:length(szB);
szB(dimB) = [];
dimsB(dimB) = [];

permuteB = permute(B, [dimsB, dimB]);
permuted_reshapedB = reshape(permuteB, prod(szB), szDimB);

tprod = permuted_reshapedA*(permuted_reshapedB.');
tprod = reshape(tprod, [szA, szB]);

end

