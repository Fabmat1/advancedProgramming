import random
import time
import numpy as np

t_start = time.time()

random.seed(143189)

Ntotal = 10000000

x, y = np.split(np.random.random(2*Ntotal), 2)

NInsideCircle = np.sum(np.square(x-0.5)+np.square(y-0.5) < 0.25)

t_end = time.time()
print("Fraction inside circle:", NInsideCircle / Ntotal)
print("Our estimation of pi", 4 * NInsideCircle / Ntotal)

print('Calculating pi took', t_end - t_start, 'seconds')
