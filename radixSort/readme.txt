Radix Sort

Sort n integers.

Perform bucket sort on every single bit(or radix number) from the lowest to the highest.

void genRandom(); generate random integers in all processes
void doRadixSort(); call bucketSort() repeatly on different bits
    void bucketSort(const DType& bv); bucket sort local data in each process
    void remap(); call Alltoall to make each process hold its corresponding elements
