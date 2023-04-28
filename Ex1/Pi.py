import random, time

t_start = time.time()

random.seed(143189)

x = []
y = []

Ntotal = 10000000

for i in range(Ntotal):
    x.append(random.random())
    y.append(random.random())

NInsideCircle = 0
for i in range(Ntotal):
    if (x[i] - 0.5) ** 2 + (y[i] - 0.5) ** 2 < 0.5 ** 2:
        NInsideCircle += 1

t_end = time.time()
print("Fraction inside circle:", NInsideCircle / Ntotal)
print("Our estimation of pi", 4 * NInsideCircle / Ntotal)

print('Calculating pi took', t_end - t_start, 'seconds')
