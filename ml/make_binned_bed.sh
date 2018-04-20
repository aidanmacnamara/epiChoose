awk ' \
    BEGIN \
    { \
        binSize = 100; \
        binIdx = 0; \
    } \
    { \
        chr = $1; \
        start = $2; \
        stop = $3; \
        for (binStart = start; binStart < (stop - binSize); binStart += binSize) { \
            print chr"\t"binStart"\t"(binStart + binSize)"\tbin-"binIdx; \
            binIdx++; \
        } \
    }' tss_win.bed \
    > tss_win_binned.bed
