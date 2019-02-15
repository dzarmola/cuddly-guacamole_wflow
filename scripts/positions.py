import json
import sys

parts = sys.argv[1].split()
outfile = sys.argv[2]

short = {
    "5a1s": [(85, 167), (304, 385)],
    "4bwz": [(35, 110), (238, 315)],
    "1zcd": [(21, 118), (229, 315)],
    "4n7w": [(10, 87), (163, 242)]
}

long = {
    "5a1s": [(63, 237), (280, 448)],
    "4bwz": [(18, 174), (224, 382)],
    "1zcd": [(9, 177), (209, 384)],
    "4n7w": [(1, 146), (152, 300)]

}

Ns = {_: [(long[_][0][0], short[_][0][0] - 1), (long[_][1][0], short[_][1][0] - 1)] for _ in long.keys()}
Cs = {_: [(short[_][0][1] + 1, long[_][0][1]), (short[_][1][1] + 1, long[_][1][1])] for _ in long.keys()}
cores = long  # short

j = {"N": Ns, "C": Cs, "core": cores}
j = {k: v for k, v in j.items() if k in parts}
with open(outfile, "w") as out:
    json.dump(j, out, sort_keys=True,
              indent=4, separators=(',', ': '))
