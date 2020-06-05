global 

Kij = bdGetValue(sys.pardef,'Kij');
aee = bdGetValue(sys.pardef,'aee');
aei = bdGetValue(sys.pardef,'aei');
aie = bdGetValue(sys.pardef,'aie');
ane = bdGetValue(sys.pardef,'ane');
ani = bdGetValue(sys.pardef,'ani');
b = bdGetValue(sys.pardef,'b');
C = bdGetValue(sys.pardef,'C');
r = bdGetValue(sys.pardef,'r');
phi = bdGetValue(sys.pardef,'phi');
Gion = bdGetValue(sys.pardef,'Gion');
Vion = bdGetValue(sys.pardef,'Vion');
thrsh = bdGetValue(sys.pardef,'thrsh');
delta = bdGetValue(sys.pardef,'delta');
VT = bdGetValue(sys.pardef,'VT');
ZT = bdGetValue(sys.pardef,'ZT');
deltaV = bdGetValue(sys.pardef,'deltaV');
deltaZ = bdGetValue(sys.pardef,'deltaZ');
I = bdGetValue(sys.pardef,'I');


[t,y] = ode23('simvec',[0:0.2:2000],[sys.vardef(1).value;sys.vardef(2).value;sys.vardef(3).value].',sys.odeoption);