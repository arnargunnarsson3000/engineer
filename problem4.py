from linear import matrix
import latex
import robotics
from sympy import symbols
import math

def insert(A,t):
    temp = A.copy()
    for i in range(A.n):
        for j in range(A.m):
            temp[i][j] = A[i][j].subs('t', t)
    return temp

def solveABC(A, b):
    q = []
    for i in range(len(b)):
        temp = [b[i], 0, 0]
        abc = A**-1 * matrix([[i] for i in temp])
        abc = [abc[0][0], abc[1][0], abc[2][0]]
        q.append(abc)

    return q

def p4():
    t = symbols('t')
    A = matrix([[t**5, t**4, t**3],[5*t**4, 4*t**3, 3*t**2],[20*t**3, 12*t**2, 6*t]])
    A0 = insert(A,0)
    print(A0)
    A2 = insert(A,2)
    print(A2)

    first = [-0.321750554396642, -4.49454236888801 + 2*math.pi, 1.61941347407016, 2.92374604209311]
    abc1 = solveABC(A2, first)
    dof = 4
    count = 1
    for i in range(len(abc1)):
        eq = "q{}{} = {}t^5 + {}t^4 + {}t^3 + {}\n".format(count, i+1, round(abc1[i][0], dof),round(abc1[i][1], dof),round(abc1[i][2], dof), round(first[i], dof))
        #eq += "dq{}{}/dt = {}t^4 + {}t^3 + {}t^2\n".format(count, i+1, round(5*abc1[i][0], dof),round(4*abc1[i][1],dof),round(3*abc1[i][2], dof))
        #eq += "d2q{}{}/dt2 = {}t^3 + {}t^2 + {}t\n".format(count, i+1, round(20*abc1[i][0], dof),round(12*abc1[i][1], dof), round(6*abc1[i][2]))
        print(eq)


    count = 2
    second = [0.26396372362570, -4.53174363933971 + 2*math.pi, 1.94807597387782, 2.96094731254482]
    abc1 = solveABC(A2, first)
    for i in range(len(abc1)):
        eq = "q{}{} = {}t^5 + {}t^4 + {}t^3 + {}\n".format(count, i+1, round(abc1[i][0], dof),round(abc1[i][1], dof),round(abc1[i][2], dof), round(first[i], dof))
        #eq += "dq{}{}/dt = {}t^4 + {}t^3 + {}t^2\n".format(count, i+1, round(5*abc1[i][0], dof),round(4*abc1[i][1],dof),round(3*abc1[i][2], dof))
        #eq += "d2q{}{}/dt2 = {}t^3 + {}t^2 + {}t\n".format(count, i+1, round(20*abc1[i][0], dof),round(12*abc1[i][1], dof), round(6*abc1[i][2]))
        print(eq)
    count = 3
    third = [-0.223476601140633, -4.55848107238232 + 2*math.pi, 2.28309001136617, 2.98768474558743]
    abc1 = solveABC(A2, second)
    for i in range(len(abc1)):
        eq = "q{}{} = {}t^5 + {}t^4 + {}t^3 + {}\n".format(count, i+1, round(abc1[i][0], dof),round(abc1[i][1], dof),round(abc1[i][2], dof), round(second[i], dof))
        #eq += "dq{}{}/dt = {}t^4 + {}t^3 + {}t^2\n".format(count, i+1, round(5*abc1[i][0], dof),round(4*abc1[i][1],dof),round(3*abc1[i][2], dof))
        #eq += "d2q{}{}/dt2 = {}t^3 + {}t^2 + {}t\n".format(count, i+1, round(20*abc1[i][0], dof),round(12*abc1[i][1], dof), round(6*abc1[i][2]))
        print(eq)
    count = 4
    fourth = [0.193621992855945, -4.57850459598843 + 2*math.pi, 2.62202212042538, 3.00770826919353]
    abc1 = solveABC(A2, third)
    for i in range(len(abc1)):
        eq = "q{}{} = {}t^5 + {}t^4 + {}t^3 + {}\n".format(count, i+1, round(abc1[i][0], dof),round(abc1[i][1], dof),round(abc1[i][2], dof), round(third[i], dof))
        #eq += "dq{}{}/dt = {}t^4 + {}t^3 + {}t^2\n".format(count, i+1, round(5*abc1[i][0], dof),round(4*abc1[i][1],dof),round(3*abc1[i][2], dof))
        #eq += "d2q{}{}/dt2 = {}t^3 + {}t^2 + {}t\n".format(count, i+1, round(20*abc1[i][0], dof),round(12*abc1[i][1], dof), round(6*abc1[i][2]))
        print(eq)
    count = 5
    fifth = [-0.170735211475283, -4.59401022420542 + 2*math.pi, 2.96352830254749, 3.02321389741052]
    abc1 = solveABC(A2, fourth)
    for i in range(len(abc1)):
        eq = "q{}{} = {}t^5 + {}t^4 + {}t^3 + {}\n".format(count, i+1, round(abc1[i][0], dof),round(abc1[i][1], dof),round(abc1[i][2], dof), round(fourth[i], dof))
        #eq += "dq{}{}/dt = {}t^4 + {}t^3 + {}t^2\n".format(count, i+1, round(5*abc1[i][0], dof),round(4*abc1[i][1],dof),round(3*abc1[i][2], dof))
        #eq += "d2q{}{}/dt2 = {}t^3 + {}t^2 + {}t\n".format(count, i+1, round(20*abc1[i][0], dof),round(12*abc1[i][1], dof), round(6*abc1[i][2]))
        print(eq)

def p5():
    t = symbols('t')
    A = matrix([[t ** 5, t ** 4, t ** 3], [5 * t ** 4, 4 * t ** 3, 3 * t ** 2], [20 * t ** 3, 12 * t ** 2, 6 * t]])
    A0 = insert(A, 0)

    A2 = insert(A, 2)

    points = [[-0.321750554396642, -4.49454236888801 + 2 * math.pi, 1.61941347407016, 2.92374604209311],
              [0.26396372362570, -4.53174363933971 + 2 * math.pi, 1.94807597387782, 2.96094731254482],
              [-0.223476601140633, -4.55848107238232 + 2 * math.pi, 2.28309001136617, 2.98768474558743],
              [0.193621992855945, -4.57850459598843 + 2 * math.pi, 2.62202212042538, 3.00770826919353],
              [-0.170735211475283, -4.59401022420542 + 2 * math.pi, 2.96352830254749, 3.02321389741052]]
    man = robotics.manipulator()
    man.addlink()
    man.addlink()
    man.addlink()
    man.addlink()
    J4 = man.H2()
    B = robotics.insert(J4, a=[0, 0, 0, 1.15], th=[None, None, math.pi, None],
               al=[-math.pi / 2, math.pi / 2, math.pi / 2, math.pi / 2], d=[1.7, 0, None, 0])
    B = [B[0][-1], B[1][-1]]
    A = matrix(5, 2)
    for i in range(len(points)):

        t1 =  B[0].subs('th1',points[i][0]).subs('th2',points[i][1]).subs('d3',points[i][2]).subs('th4',points[i][3])
        t2 =  B[1].subs('th1',points[i][0]).subs('th2',points[i][1]).subs('d3',points[i][2]).subs('th4',points[i][3])

        A[i][0] = round(t1, 3)
        A[i][1] = round(t2, 3)
    ini = [i for i in A[:-1]]
    fin = [i for i in A[1:]]
    r = 3
    for i in range(len(fin)):
        q = solveABC(A2, fin[i])
        str = "p_x(t) = {}t^5 + {}t^4 + {}t^3 + {}\\\\\n".format(round(q[0][0],r), round(q[0][1],r), round(q[0][2],r),ini[i][0])
        str += "p_y(t) = {}t^5 + {}t^4 + {}t^3 + {}\\\\\n".format(round(q[1][0],r), round(q[1][1],r), round(q[1][2],r),ini[i][1])
        print(str)

def p6():
    man = robotics.manipulator()
    man.addlink('r')
    man.addlink('r')
    man.addlink('p')
    man.addlink('r')
    J = man.Jacobian()
    print(J)
    B = robotics.insert(J, a=[0, 0, 0, 1.15], th=[None, None, math.pi, None],
                        al=[-math.pi / 2, math.pi / 2, math.pi / 2, math.pi / 2], d=[1.7, 0, None, 0])

    print(latex.Lmatrix(B.simplify(), num=['al','a','th','d']))



if __name__ == "__main__":
    #p4()
    #p5()
    p6()