import time

import numpy

numpy.random.seed(143189)

Ntotal = 10000000

t0 = time.time()

ArrayOfNumbers = numpy.random.random(Ntotal)

Min = numpy.min(ArrayOfNumbers)
Max = numpy.max(ArrayOfNumbers)
Mean = numpy.mean(ArrayOfNumbers)

t1 = time.time()

print("The smallest number is",Min)
print("The largest number is",Max)
print("The mean value is",Mean)

print(f"Time to run: {round(t1-t0, 3)}")
