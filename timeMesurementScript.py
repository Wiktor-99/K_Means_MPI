import subprocess
import datetime

retries = 10
processes = 20

iterations = [10, 50, 100, 250, 500, 750, 1000]
centroid = [3, 4, 5, 6, 7, 8, 9, 10]

allTimes = []
for iteration in iterations:
    times = []
    for numberOfProcesses in range(processes):
        elapsed_time = 0
        for i in range(retries):
            st = datetime.datetime.now()
            subprocess.run(['mpiexec', '-n', str(numberOfProcesses), 'a.out', '7', str(iteration)])
            et = datetime.datetime.now()
            elapsed_time += (et - st).seconds
        times.append(elapsed_time / retries)
    allTimes.append(times)

print('Times for different iterations and processes')
print(allTimes)

allTimes = []
for c in centroid:
    times = []
    for numberOfProcesses in range(processes):
        elapsed_time = 0
        for i in range(retries):
            st = datetime.datetime.now()
            subprocess.run(['mpiexec', '-n', str(numberOfProcesses), 'a.out', str(c), '100'])
            et = datetime.datetime.now()
            elapsed_time += (et - st).seconds
        times.append(elapsed_time / retries)
    allTimes.append(times)

print('Times for different centroids and processes')
print(allTimes)
