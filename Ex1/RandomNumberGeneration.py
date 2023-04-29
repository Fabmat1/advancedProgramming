import random
import time
random.seed(143189)
ListOfNumbers = []
Ntotal = 10000000

print("We are now drawing",Ntotal," random numbers between 0.0 and 1.0 from a uniform distribution.")

t0 = time.time()

for i in range(Ntotal):
    ListOfNumbers.append(random.random())


Max = ListOfNumbers[0]
Min = ListOfNumbers[0]

Mean = 0.0

for i in range(Ntotal):
    if ListOfNumbers[i]>Max:
        Max = ListOfNumbers[i]

for i in range(Ntotal):
    if ListOfNumbers[i]<Min:
        Min = ListOfNumbers[i]

for i in range(Ntotal):
    Mean = Mean + ListOfNumbers[i]/Ntotal


t1 = time.time()

print("The smallest number is",Min)
print("The largest number is",Max)
print("The mean value is",Mean)

print("The smallest number is (calculated by the min function)",min(ListOfNumbers))
print("The largest number is(calculated by the max function)",max(ListOfNumbers))
print(f"It took {round(t1-t0, 3)} seconds")
