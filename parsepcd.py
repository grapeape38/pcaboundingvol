import sys
out = open("testpointdata.js", "w+")
for i in range(1, len(sys.argv)):
    fname = sys.argv[i]
    objname = fname.split('.')[0]
    f = open(fname, "r")
    fi = f.readlines()

    lines = fi[:]
    lc = 0
    for line in lines:
        lc+=1

    skip = (lc-10) / 2000

    lines = fi[:]
    out.write("var " + objname + " = [")
    i = 0
    for line in lines:
        if i > 10 and i % skip == 0:
            coords = line.split(' ')
            out.write("[{},{},{}],\n".format(coords[0], coords[2].strip(), coords[1]))
        i+=1
    out.write('];\n')
    f.close()

out.close()


