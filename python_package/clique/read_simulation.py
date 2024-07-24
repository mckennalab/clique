import gzip

def read_n(filename, x):
    with gzip.open(filename, 'rt') as fh:
        while True:
            data = [fh.readline().strip() for _ in range(x)]

            if len(data[0]) == 0:
                break

            yield data


output = open("read_assignments.txt", "w")
output.write("name\ttag\tchimera\tlength\terror_free\tidentity\n")
for count, nlines in enumerate(read_n("reads.fastq.gz", 4)):
    # print(nlines[0])
    t = nlines[0].split(" ")
    length = 0
    error_free_len = 0
    identity = 0.0
    chimera = False
    if "chimera" in nlines[0]:
        chimera = True
    for tag in t:
        if tag.startswith("length="):
            length = int(tag.split("=")[1])
        if tag.startswith("error"):
            error_free_len = int(tag.split("=")[1])
        if tag.startswith("read_identity"):
            identity = float(tag.split("=")[1][0:len(tag.split("=")[1]) - 1])
    else:
        output.write(
            t[0].lstrip("@") + "\t" + t[1].split(",")[0] + "\t" + str(chimera) + "\t" + str(length) + "\t" + str(
                error_free_len) + "\t" + str(identity) + "\n")
output.close()


