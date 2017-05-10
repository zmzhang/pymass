cuda_add_library(gsort STATIC cuda_sort.h cuda_sort.cu)
cuda_add_executable(cuda_sort cuda_sort.cu test.cu)