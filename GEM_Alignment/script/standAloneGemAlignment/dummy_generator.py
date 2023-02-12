

f = open("ME11_misalignment_dx.csv", "w")

x=1
y=0
z=0
phix=0
phiy=0
phiz=0

for i in [-1, 1]:
  for j in range (1,37):
    if j <10:
      f.writelines("{i}0{j}, {x}, {y}, {z}, {phix}, {phiy}, {phiz}\n".format(i=i, j=j, x=x, y=y, z=z, phix=phix, phiy=phiy, phiz=phiz))
    else:
      f.writelines("{i}{j}, {x}, {y}, {z}, {phix}, {phiy}, {phiz}\n".format(i=i, j=j, x=x, y=y, z=z, phix=phix, phiy=phiy, phiz=phiz))

f.close()
