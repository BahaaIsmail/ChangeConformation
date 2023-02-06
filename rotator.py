
from math import pi , sqrt , cos , sin , atan2
import sys
### given the following 
# the axese of rotation in terms of two points p1 and p2 
# the point p which will be rotated
# the angle t through whih the point will be rotated 

### notes
# the direction of the rotation vector and the direction of rotation are carefully adjusted to give the right values 






#__________________________________________________________________________
#  extracting atom numbers , connectivities and cartesian coordinates from the original pdb file 
#  the entire "ATOM" section will be restorerd to be exported to the new pdb file
pdb = 'seq.pdb'

sf = open(pdb)

atomsymb = {}
cc = {}
conect = {}
ATOM = {}
k = 1
for line in sf : 
    if 'ATOM' in line : 
        line = line.split()
        ATOM[k] = line
        xyz = line[5:8]
        cc[k] = [float(i) for i in xyz]
        atomsymb[k] = line[-1] 
        k += 1
    elif 'CONECT' in line : 
        line = line.split()
        if len(line) > 2 : 
            atoms = [int(i) for i in line[1:]]
            conect[atoms[0]] = atoms[1:]     
sf.close()

#__________________________________________________________________________
# extracting the torsional angles to be rotated 
sf = open('rotinp')
torsionals = []
for line in sf :   
    line = line.split()
    if len(line) == 5 : 
      torsionals += [[int(i) for i in line[:4] + [float(line[4])]]]

sf.close()



#_______________________________________________________________
### functions to calculate the values of the torsional angles  
def dot(v1,v2) : 
    return sum([v1[i]*v2[i] for i in range(len(v1))])

def norm (v) : 
    return sqrt(dot(v,v))  

def disp(a1,a2): 
    return [a2[i] - a1[i] for i in range(len(a1))]    
    
def dist(a1,a2): 
    v = disp(a1,a2)
    return norm(v)

def uv(v) : 
    n = norm(v)
    return [i/n for i in v]

def cross (v1,v2) :  # v1 x v2
    return [v1[1]*v2[2]-v1[2]*v2[1] , v1[2]*v2[0]-v1[0]*v2[2] , v1[0]*v2[1]-v1[1]*v2[0]]    

print ()

#_______________________________________________________________
### calculating the values of the current torsional angles and the the ritation value  
radeg = 180/pi   # convet from radians to degrees
for t in range(len(torsionals)) : 
    [i,j,k,l] = torsionals[t][:4]
    vij , vjk , vkl = uv(disp(cc[j],cc[i])) , uv(disp(cc[k],cc[j])) , uv(disp(cc[l],cc[k]))
    ni = cross(vij,vjk)
    nk = cross(vjk,vkl)
    m  = cross(ni,vjk)
    x , y = dot(ni,nk) , dot(m ,nk)
    value = atan2(y,x)*radeg
    print ('the angle  ' , i , j , k , l , ' =  ' , round(value,1) , ' ---> '  , torsionals[t][4] ,' <---- ' , round(torsionals[t][4]-value,1))
    value = torsionals[t][4]-value
    torsionals[t] += [round(value,2)]
print()

#__________________________________________________________________________
# finding out the group of atoms associated with each rotatble torsional angle 
# a list of the torsional angles to be rotated with the desired angle 
#    is provided in the inpit file "rotinp"
def diff (L1 , L2) : 
    return [i for i in L1 if not i in L2]    

def tribe(i,j)  :     
    itribe = [k for k in conect[i] if not k == j]
    jtribe = [k for k in conect[j] if not k == i]
    updated = 1
    while updated : 
        updated = 0   
        iupdate = itribe + []
        for k in iupdate : 
            if  not  k in jtribe : 
                update = diff(conect[k],itribe+[i,j]) 
                if update :  
                        itribe += update
                        updated = 1

        jupdate = jtribe + []
        for k in jupdate : 
            if  not  k in itribe : 
                update = diff(conect[k],jtribe+[i,j]) 
                if update :     
                        jtribe += update
                        updated = 1

    shared = [k for k in itribe if k in jtribe]    
    itribe = diff(itribe,shared)
    jtribe = diff(jtribe,shared)   
    ni = len(itribe)
    nj = len(jtribe)
    if ni <= nj : 
        return -1 , sorted(itribe)
    else : 
        return  1 , sorted(jtribe)


### loopiong over the torsionals to determine the atomgrups to be rotated
atomgroups = []
for t in range(len(torsionals)) :
    [i,j] = torsionals[t][1:3]
    d , atomgroup = tribe(i,j)
    torsionals[t][-1] = d*torsionals[t][-1]
    atomgroups += [atomgroup]



### the rotating function 
def prepare(p1,p2,t) : 
    
    [x , y , z] = p1

    T   = [[ 1 , 0 , 0 , -x ] , [ 0 , 1 , 0 , -y ] , [ 0 , 0 , 1 , -z ] , [ 0 , 0 , 0 ,  1 ] ] 
    T_1 = [[ 1 , 0 , 0 ,  x ] , [ 0 , 1 , 0 ,  y ] , [ 0 , 0 , 1 ,  z ] , [ 0 , 0 , 0 ,  1 ] ] 

    v =  [p2[i]-p1[i] for i in range(3)] 
    n = sqrt(sum([i**2 for i in v]))
    u = [i/n for i in v]
    [a , b , c] = u  
    d = sqrt(b**2+c**2)   
    b , c = b/d , c/d ; 

    Rx   = [[ 1 , 0 ,  0 , 0 ] , [ 0 , c , -b , 0 ] , [  0 ,  b , c , 0 ] , [ 0 , 0 , 0 , 1 ] ] 
    Rx_1 = [[ 1 , 0 ,  0 , 0 ] , [ 0 , c ,  b , 0 ] , [  0 , -b , c , 0 ] , [ 0 , 0 , 0 , 1 ] ] 
    Ry   = [[ d , 0 , -a , 0 ] , [ 0 , 1 ,  0 , 0 ] , [  a ,  0 , d , 0 ] , [ 0 , 0 , 0 , 1 ] ]  
    Ry_1 = [[ d , 0 ,  a , 0 ] , [ 0 , 1 ,  0 , 0 ] , [ -a ,  0 , d , 0 ] , [ 0 , 0 , 0 , 1 ] ] 

    c , s = cos(t) , sin(t)
    Rz = [[ c , -s , 0 , 0] , [ s ,  c , 0 , 0] , [ 0 ,  0 , 1 , 0] , [ 0 ,  0 , 0 , 1] ] 


    def MM (M1,M2): 
        m , n = len(M1) , len(M2)
        M = [[0 for j in range(n)] for i in range(m)]
        for i in range(m) : 
          for j in range(n) : 
            M[i][j] = sum([M1[i][k]*M2[k][j] for k in range(m)])
        return M
     

    R = MM(Rx,T)
    R = MM(Ry,R)
    R = MM(Rz,R)
    R = MM(Ry_1,R)
    R = MM(Rx_1,R)
    R = MM(T_1,R)
    return R


def rotate(p,R):
    p += [1]
    f = [[p[i]] for i in range(len(p))]

    p = []
    for i in range(4) : 
      for j in range(len(f[0])) : 
        p += [sum([R[i][k]*f[k][j] for k in range(4)])]
    #print p[:3] , '________________'
    return p[:3]


### looping over the torsionals to be rotated 
n = len(torsionals)
for k in range(n) :        
    [i,j] = torsionals[k][1:3]
    p1 , p2 = cc[i] , cc[j]    
    t = torsionals[k][-1]*pi/180
    atomgroup = atomgroups[k]
    
    R = prepare(p1,p2,t)
    for l in atomgroup : 
        p = cc[l]+[]
        p = rotate(p,R) 
        cc[l] = p + []





### exporting the pdb file
tf = open('r_'+pdb,'w')
sf = open(pdb)



for i in ATOM : 
    ATOM[i][5:8] = [str(round(j,3)) for j in cc[i]]

k = 0
for line in sf : 
    if 'ATOM' in line :
        k += 1
        L = ATOM[k]
        tf.write('%4s%7s%5s%2s%8s%12s%8s%8s%6s%6s%12s\n' %(L[0],L[1],L[2],L[3],L[4],L[5],L[6],L[7],L[8],L[9],L[10]))
    else : 
      tf.write(line)


sf.close()
tf.close()
