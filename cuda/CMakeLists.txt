cuda_add_library(gsort STATIC cuda_sort.h cuda_sort.cu
                 OPTIONS -DPYMASS_STATIC_DEFINE)
cuda_add_executable(cuda_sort cuda_sort.cu test.cu
                 OPTIONS -DPYMASS_STATIC_DEFINE)