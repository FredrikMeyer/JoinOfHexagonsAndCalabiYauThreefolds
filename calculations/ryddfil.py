

with open("alltriangs.txt") as f:
    L = list(f.readline())

n = len(L)
for i in range(n):
    if L[i] in ["[","(", "<"]:
        L[i] = "{"
    elif L[i] in ["]", ")", ">"]:
        L[i] = "}"

S = "".join(L)

with open("alltriangsm2.txt", "w") as f:
    f.write(S)
