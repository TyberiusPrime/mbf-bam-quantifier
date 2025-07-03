a = [x.strip() for x in open("../input_BD_CLS1.txt", "r").readlines() if x.strip()]
b = [x.strip() for x in open("../input_BD_CLS2.txt", "r").readlines() if x.strip()]
c = [x.strip() for x in open("../input_BD_CLS3.txt", "r").readlines() if x.strip()]

should = []
for ax in a:
    for bx in b:
        for cx in c:
            should.append(f"{ax}_{bx}_{cx}")
with open("valid.barcode",'w') as op:
    for x in should:
        op.write(f"{x}\n")
