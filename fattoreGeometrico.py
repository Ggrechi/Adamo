import multiprocessing
from multiprocessing.pool import Pool
import time

cpu = multiprocessing.cpu_count() 
#print cpu
a = float(raw_input("dimensione scintillatore 1 (m) a: "))
b = float(raw_input("dimensione scintillatore 2 (m) b: "))
d = float(raw_input("distanza scintillatori (m): "))

def partialSum(v):
    [N_min,N_max,N,a,b,d] = v
    #print range(N_min,N_max)
    pSum = float(0.0)
    for i in range(N_min,N_max):
        #print sum
        for j in range(N):
            for n in range(N):
                for m in range(N):
                    pSum += 1/(d**2 + (a*n/N-b*j/N)**2 + (a*m/N-b*i/N)**2)**2
    return pSum


#for N in range(1,int(raw_input("numero massimo di suddivisioni finale: "))+1):
"""
#if __name__ == '__main__':
    #N=int(raw_input("numero di intervalli N: "))
"""
    startTime = time.clock()
    p = Pool(cpu)
    result = float(0.0)
    factor = ((a*b*d/(N**2))**2)
    v = [[N*spac/cpu,N*(spac+1)/cpu,N,a,b,d] for spac in range(cpu)]

    res = p.map(partialSum, v)

    p.close()
    p.join()
    #for index in range(len(res)):
    result += factor * sum(res)
    elapsedTime = time.clock() - startTime

    print N,result,result*80,"Hz\t",elapsedTime
