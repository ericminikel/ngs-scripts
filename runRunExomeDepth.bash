# bams: /humgen/atgu1/fs03/eminikel/sandbox/bamlist.txt
# outdir: /humgen/atgu1/fs03/eminikel/sandbox/output

# submit with:
# bsub -o /humgen/atgu1/fs03/eminikel/sandbox/runRunExomeDepthBsubOutput.out \
#     bash /humgen/atgu1/fs03/eminikel/sandbox/runRunExomeDepth.bash

# runExomeDepth.r is in my src/ directory which I have in path

runExomeDepth.r -b /humgen/atgu1/fs03/eminikel/sandbox/bamlist.txt \
                  -o /humgen/atgu1/fs03/eminikel/sandbox/output \
                  -v > /humgen/atgu1/fs03/eminikel/sandbox/runExomeDepthOutput.txt
