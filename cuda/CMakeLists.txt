cuda_add_library(gsort STATIC cuda_sort.h cuda_sort.cu
                 OPTIONS -DPYMASS_STATIC_DEFINE)