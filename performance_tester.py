import os
import math as m
import subprocess
import resource
import os
import random

SEED = 1.3

def main() :
    distances = [0.1, 0.05, 0.01, 0.005, 0.001]
    points = []
    s = SEED
    points = [int(s)]
    for _ in range(49) :
        s = SEED*s
        points.append(int(s))
    test_unif(points, distances)

def test_unif(points, distances) :
    with open("perf_graphes/unif_wsl.csv", 'w') as f :
        f.write("d, n, t, m\n")
        for i in range(len(distances)) :
            res = []
            for j in range(len(points)) :
                d = distances[i]
                n = points[j]
                print("test avec", n, "points et une distance", d)
                points_generator_unif(n, d)
                p = subprocess.run(["/usr/bin/time", "-f%e", "python3", "main.py", "test_points.pts"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                t = float(p.stderr)
                p = subprocess.run(["/usr/bin/time", "-f%M", "python3", "main.py", "test_points.pts"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                m = float(p.stderr)
                f.write(str(d) + ", " + str(n) + ", " + str(t) + ", " + str(m) + "\n")
                print(str(i*len(points) + (j+1)) + " / " + str(len(distances)*len(points)) + " tests completed")
    os.remove("test_points.pts")

def points_generator_unif(n, d) :
    with open("test_points.pts", 'w') as file :
        file.write(str(d) + "\n")
        for i in range(n) :
            file.write(str(random.random()) + ", " + str(random.random()) + '\n')

if __name__ == '__main__':
    main()
