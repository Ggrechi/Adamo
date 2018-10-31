from scipy.integrate import nquad
#d=distanza tra i rivelatori
#x1,y1 lati x y del rivelatore 1
#x2,y2 lati x y del rivelatore 2
d=0.04
x_1=0.137
y_1=0.115
x_2=0.137
y_2=0.115

def function(x,y,a,b):
	return d**2/(d**2+(a-x)**2+(b-y)**2)**2

result, result_err = nquad(function,[[0,x_1], [0,y_1], [0,x_2], [0,y_2]])

print result
print result_err
print result*80,"Hz"
