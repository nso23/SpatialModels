#!/usr/bin/env python3

import subprocess
import sys
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import curve_fit
import sympy
import random

def myRun(runFile, target):
    run = subprocess.Popen(["/usr/local/bin/slim", runFile],
                        stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                        universal_newlines = True)
    stdout = []
    results = []
    count = 0
    res = [[0.0, 0.0], 0.0]
    while True:
        count += 1
        line = run.stdout.readline()
        if "#N:" in line:
            results.append(int(line.split("#N: ")[1]))
            if results[-1] > target*1.5:
                return [0,results[-1]]
                
        if count > 100:
            if results[-1] == 0:
                break
            try:
                res = fitCurve(results)
                if res[1] > 0.9:
                    if abs(np.mean(results[-10:-1]) - res[0][1])/res[0][1] < 0.01:
                        break
                if resuls[-1] == 0:
                    break
            except:
                continue
        if count > 5000:
            res[0] = 0
            break
    run.kill()
    return res[0]

def func(x, a, k, q, b):
        return a + k-a/((1+q*scipy.special.expit(-b*x)))

def fitCurve(results):
        Nr = np.array([sympy.N(i) for i in results],dtype='float64')
        params, params_covariance = curve_fit(func, [i for i in range(0, len(Nr))], Nr,p0=[0, 10000, 1, 1])
        N_fit = []
        for i in [i for i in range(0, len(Nr))]:
            N_fit.append(func(i, params[0], params[1], params[2], params[3]))
        ss_res = np.sum((Nr - N_fit) ** 2)
        ss_tot = np.sum((Nr - np.mean(Nr)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        return [params, r2]

def minFunc(outFile, target):
    params = myRun(outFile, target)
    error = (params[1]-target)/target
    return min(error, 1)

def chooseM(score):
    if len([i for i in score[1] if i < 1]) < 3:
        return [random.uniform(0, 0.1),random.uniform(0, 0.1)]
    y1, x1, z1 = (list(t) for t in zip(*sorted(zip(score[1], score[0], score[2]))))
    x2, y2, z2 = (list(t) for t in zip(*sorted(zip(score[0], score[1], score[2]))))
    i = x2.index(x1[0])
    minv = x2[i]
    maxv = x2[i]
    if z1[0] > 0:
            if z2[i] - z2[i-1] < 0:
                return [x2[i]/2 + x2[i+1]/2]
            else:
                if i > 0:
                    return [x2[i]/2 + x2[i-1]/2]
                else:
                    return [x2[i]*3/4]
    else:
        if z2[i] - z2[i-1] > 0:
            return [x2[i]/2 + x2[i+1]/2]
        else:
            return [x2[i]/2 + x2[i-1]/2]

def createSLiMFile(inFile, outFile, m):
    newFile = []
    target = 1
    with open(inFile, "r") as f:
        for line in f:
            if "defineConstant(\"m\"" in line:
                newFile.append("defineConstant(\"m\", " + str(m) + ");\n")
            else:
                newFile.append(line)
            if "defineConstant(\"A\"" in line:
                tmp = line.split(",")[1]
                target *= int(tmp.split(")")[0])
            if "defineConstant(\"Rd\"" in line:
                tmp = line.split(",")[1]
                target *= int(tmp.split(")")[0])
    with open(outFile, "w+") as f:
        for line in newFile:
            f.write(line)
    return target

def findM(args):
    inFile = args[1]
    outFile = "outfile.slim"
    score = [[],[],[]]
    error = 1.0
    while error > 0.01:
        m = chooseM(score)
        target = createSLiMFile(inFile, outFile, m[0])
        res = minFunc(outFile, target)
        score[0].append(m[0])
        score[1].append(abs(res))
        score[2].append(res)
        res = abs(res)
        if res < error:
            error = res
        print(str(m[0]) + " " + str(res))
    y1, x1= (list(t) for t in zip(*sorted(zip(score[1], score[0]))))
    return x1[0]

if __name__ == "__main__":
    # Get the input file containing the parameoutFileters
    m = findM(sys.argv)
    createSLiMFile(sys.argv[1], "outfile.slim", m)
