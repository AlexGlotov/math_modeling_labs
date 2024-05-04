from math import exp, pi, cos
import numpy as np

t = []
points = []


with open('data7.txt', 'r') as f:
    data = f.readlines()

tau, step, N = data[0].split()
tau = float(tau)
step = float(step)
N = int(N)


for row in data[1:]:
    T, *U = row.split()
    U = [float(el) for el in U]
    points.append(U)
    t.append(float(T))
    
points = np.array(points)
real_points = np.zeros( (len(t) , N + 1))

k = 0
for time in t:
    x = 0.0
    for j in range(N + 1):  
        res = x*time*time
        for i in range(0, 50):
            l = pi * (1 / 2 + float(i))
            res += 2 / (l**6) * (exp(-l*l*time) + l*l*time - 1) * cos(l*x)
        real_points[k,j] = res
        x += step
    k += 1

print(np.max(abs(points - real_points)))
