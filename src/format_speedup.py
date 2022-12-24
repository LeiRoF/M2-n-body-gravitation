import regex as re

times = []
with open('data/tmp_speedup.txt', 'r') as f:
    for line in f:
        res = re.search(r'real\s([0-9]*)m([0-9]*\.[0-9]*)s', line)

        if res is None:
            continue

        time = float(res.group(1)) * 60 + float(res.group(2))
        times.append(str(time))

with open('data/speedup.txt', 'a') as f:
    f.write("\t".join(times) + "\n")