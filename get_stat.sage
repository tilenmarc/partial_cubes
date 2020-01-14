from sage.all_cmdline import *

load('PCproperties.sage')

o = open('partial_cubes_i5.txt')

results = {}
for i in range(1, 7):
    results[i] = 14 * [0]
    # all, well-embedded, com, om, pasch, hypercellular, Polat, cellular, median, peano, netlike,  treezone, LOP, almost median

index = 0
for line in o:
    index += 1
    # if index % 100 == 0:
    #     print index
    G = Graph(line)
    if G.order() == 0:
        continue
    PC = PCproperties(G)

    results[PC.dimension][0] += 1
    if PC.is_well_embedded():
        results[PC.dimension][1] += 1
        if PC.is_COM():
            results[PC.dimension][2] += 1
            if PC.is_OM():
                results[PC.dimension][3] += 1
            if PC.is_pasch():
                results[PC.dimension][4] += 1
                if PC.is_hypercellular():
                    results[PC.dimension][5] += 1
                    if PC.is_polat():
                        results[PC.dimension][6] += 1
                        if PC.is_cellular():
                            results[PC.dimension][7] += 1
                        if PC.is_median():
                            results[PC.dimension][8] += 1
                if PC.is_peano():
                    results[PC.dimension][9] += 1
                if PC.is_tree_zone():
                    results[PC.dimension][11] += 1
                if PC.is_netlike():
                    results[PC.dimension][10] += 1
            if PC.is_almost_median():
                results[PC.dimension][13] += 1
                if PC.is_LOP():
                    results[PC.dimension][12] += 1
o.close()

for i in range(1, 7):
    print "isometric dimension " + str(i)
    print("all " + str(results[i][0]) + ", well-embedded " + str(results[i][1]) + ", com " + str(results[i][2]) + ", om " + str(results[i][3]) + ", pasch " + str(results[i][4]) + ", hypercellular " + str(results[i][5]) + ", Polat " + str(results[i][6]) + ", cellular " + str(results[i][7]) + ", median " + str(results[i][8]) + ", peano " + str(results[i][9]) + ", netlike " + str(results[i][10]) + ", treezone " + str(results[i][11]) + ", LOP " + str(results[i][12]) + ", almost median " + str(results[i][13]))
