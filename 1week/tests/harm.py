def f(i):
    return 1.0/i

n = 10000000
s = 0
for i in range(1,n+1):
    s += f(i)
print("s=",s)

