#------------------------------ Small Useful Functions ------------------------------

def find_closest_root(roots, a_root):
    '''Among roots, returns the closest one to a_root'''
    r = roots[0]    
    ab = abs(r-a_root)
    for root in roots:
        nab = abs(root-a_root)
        if nab < ab:
            r = root
            ab = nab
    return r

# 3 digits precision, only for output
n3 = lambda x: n(x,digits=3)
n6 = lambda x: n(x,digits=6)

# is_infinity for RIF
is_infinity = lambda x: RIF(x).upper()==+infinity or RIF(x).lower()==-+infinity

# radius value by name
r_value = lambda x:I if x=='1' else r if x=='r' else s

#radius name by value
r_name = lambda x:'1' if x>r else 'r' if x>s else 's'

# mq value by q
mq_value = lambda x:m1 if x>r else mr if x>s else ms

# zq value by q
zq_value = lambda x:z1 if x>r else zr if x>s else zs

# triangle type detection
triangle_name = lambda t: r_name(t[0])+r_name(t[1])+r_name(t[2])
# triangle type detection, sorted
sorted_triangle_name = lambda t: ''.join(sorted(triangle_name(t)))


# triangle type detection
triangle_value = lambda t: (r_value(t[0]),r_value(t[1]),r_value(t[2]))


# area of a triangle with sides of length a, b and c (all of them are intervals or Q)
area = lambda a,b,c: sqrt((a+b+c)*(a+b-c)*(a-b+c)*(b+c-a))/QQ(4)

# area of a tight triangle with discs of radii ra, rb and rc (all of them are intervals or Q)
d_area = lambda ra,rb,rc: area(ra+rb,rb+rc,ra+rc)

# angle opposite to edge a, (all of them are intervals or Q)
angle = lambda a,b,c: arccos(RIF((b^2+c^2-a^2)/(QQ(2)*b*c)))

# symbolic,no RIF angle opposite to edge a, (all of them are intervals or Q)
symb_angle = lambda a,b,c: arccos((b^2+c^2-a^2)/(QQ(2)*b*c))


# given t=(ra,rb,rc) returns the angle in the center of the circle rb (all of them are intervals or Q)
tight_angle = lambda t: angle(t[0]+t[2], t[0]+t[1], t[1]+t[2])

# L(ri,rj,rk): lower bound on angle in rj in any triangle which can appears in a FM-triangulation of a saturated packing (ri, rj and rk: sizes of circles)
lower_bound_angle = lambda ri,rj,rk: RIF(arcsin(RIF(min(ri/(rj+ri+2*s),rk/(rj+rk+2*s)))))

# area of the triangle with sides of length a, b and c covered by the discs of radius ra, rb and rc
cov = lambda a,b,c,ra,rb,rc: angle(a,b,c)*ra^2/QQ(2)+angle(b,c,a)*rb^2/QQ(2)+angle(c,a,b)*rc^2/QQ(2)
symb_cov = lambda a,b,c,ra,rb,rc: symb_angle(a,b,c)*ra^2/QQ(2)+symb_angle(b,c,a)*rb^2/QQ(2)+symb_angle(c,a,b)*rc^2/QQ(2)


# excess of the same triangle
E=lambda a,b,c,ra,rb,rc: d_opt*area(a,b,c)-cov(a,b,c,ra,rb,rc)
symb_E = lambda a,b,c,ra,rb,rc: d_opt*area(a,b,c)-symb_cov(a,b,c,ra,rb,rc)

# given q, returns all triples of radii with q in the middle
tripl = lambda q: [(a,q,b) for (a,b) in [(I,I),(I,r),(I,s),(r,r),(r,s),(s,s)]]


# value a tq f(a)f(a+eps)<0 (point of the sign changing if exists)
def sign_change_point(f,a,b,EPS):
    if (b-a)<EPS:
        return a
    c=(a+b)/QQ(2)
    if f(a)*f(c)<QQ(0):
        return sign_change_point(f,a,c,EPS)
    else:
        return sign_change_point(f,c,b,EPS)

#---------------------------- Radius of the support circle
# 2% faster ;-)

# TODO: problem radius is infinity since A nd B contain 0
# return (nb quasi-tight, nb not feasible, nb ok)

Ss = lambda a,b,c: (a+b+c)*(a+b-c)*(a-b+c)*(b+c-a)
Aa = lambda a,b,c,ra,rb,rc: RIF(a^4 - 2*a^2*b^2 + b^4 + c^4 + 4*a^2*ra^2 + 4*b^2*rb^2 + 4*c^2*rc^2 - 2*(a^2 + b^2)*c^2 - 4*(a^2 + b^2 - c^2)*ra*rb - 4*((a^2 - b^2 + c^2)*ra - (a^2 - b^2 - c^2)*rb)*rc)
Bb = lambda a,b,c,ra,rb,rc: 4*a^2*ra^3 + 4*b^2*rb^3 + 4*c^2*rc^3 - 2*(a^2 + b^2 - c^2)*ra*rb^2 - 2*((a^2 - b^2 + c^2)*ra - (a^2 - b^2 - c^2)*rb)*rc^2 + 2*(a^4 - a^2*b^2 - a^2*c^2)*ra - 2*(a^2*b^2 - b^4 + b^2*c^2 + (a^2 + b^2 - c^2)*ra^2)*rb + 2*(c^4 - (a^2 + b^2)*c^2 - (a^2 - b^2 + c^2)*ra^2 + (a^2 - b^2 - c^2)*rb^2)*rc
Cc = lambda a,b,c,ra,rb,rc: a^2*b^2*c^2 + a^2*ra^4 + b^2*rb^4 + c^2*rc^4 + (a^4 - a^2*b^2 - a^2*c^2)*ra^2 - (a^2*b^2 - b^4 + b^2*c^2 + (a^2 + b^2 - c^2)*ra^2)*rb^2 + (c^4 - (a^2 + b^2)*c^2 - (a^2 - b^2 + c^2)*ra^2 + (a^2 - b^2 - c^2)*rb^2)*rc^2
Dd = lambda a,b,c,ra,rb,rc: 4*(Ss(a,b,c))*(a+rb-rc)*(a-rb+rc)*(b+ra-rc)*(b-ra+rc)*(c+ra-rb)*(c-ra+rb)

def radius(a,b,c,ra,rb,rc):
    S=Ss(a,b,c)
    A=Aa(a,b,c,ra,rb,rc)
    B=Bb(a,b,c,ra,rb,rc)
    C=Cc(a,b,c,ra,rb,rc)
    D=Dd(a,b,c,ra,rb,rc)
    return radiusABCD(a,b,c,ra,rb,rc,A,B,C,D)

def radiusABCD(a,b,c,ra,rb,rc,A,B,C,D):
    S=Ss(a,b,c)
    if A.contains_zero():
        if B.contains_zero() or D.contains_zero():
            raise NameError("A and (B or C) contain zero in radiusABCD")
        return -C/B-A*C^2/sqrt(D)^3
    else:
        r1=(-B-sqrt(D))/(QQ(2)*A)
        r2=(-B+sqrt(D))/(QQ(2)*A)
        return r1 if ((r1>0 and r2>r1) or (r2<0)) else r2
  
def dist(a,b,c,ra,rb,rc,R):
    ''' Signed distance from the center of the support circle of radius R to the edge a '''
    S=Ss(a,b,c)
    return ((b^2+c^2-a^2 + 2*R*(rb+rc-2*ra)+rb^2+rc^2-2*ra^2)*a^2 + (2*R+rb+rc)*(rb-rc)*(b^2-c^2))/(2*a*sqrt(S))

#--------------------------------------- Coronas generator ---------------

# check whether the angle vector k does indeed correspond to a corona
def coronable(k):
    # check that there is an integer number of discs in the corona
    a=mod(k[1]+k[2],2)==0 and mod(k[1]+k[4],2)==0 and mod(k[2]+k[4],2)==0
    # if yxy, then either only y's or yxz or yxx
    b=k[1]!=0 or k[2]!=0 or k[0]==0 or (k[1]==0 and k[2]==0 and k[3]==0 and k[4]==0 and k[5]==0)
    c=k[1]!=0 or k[4]!=0 or k[3]==0 or (k[0]==0 and k[1]==0 and k[2]==0 and k[4]==0 and k[5]==0)
    d=k[2]!=0 or k[4]!=0 or k[5]==0 or (k[0]==0 and k[1]==0 and k[2]==0 and k[3]==0 and k[4]==0)
    # cycles on the triangle between 1,r,s
    z=min(k[1],k[2],k[4])
    e=mod(k[1]-z,2)==0 and mod(k[2]-z,2)==0 and mod(k[4]-z,2)==0
    return a and b and c and d and e

# lexicographic enumeration of vectors x such that, weighted by u, the sum is at most U
def odo(u,U):
    x=[0 for i in u]
    z=0
    fini=false
    while not fini:
        yield x
        i=0
        while i<len(u) and z+u[i]>U:
            z-=x[i]*u[i]
            x[i]=0
            i+=1
        if i<len(u):
            x[i]+=1
            z+=u[i]
        else:
            fini=true
            
def corona_generator(q):
    ttriangles_q = tripl(q)
    min_angles_q = [lower_bound_angle(*t) for t in ttriangles_q] # pas tight !!!
    for vec in odo(min_angles_q, 2*RIFpi):
        if coronable(vec):
            res = []
            for (triangle, multiplicity) in zip(ttriangles_q, vec):
                res+=[triangle]*multiplicity
            yield (res)

def false_corona_generator(q):
    # using a greater value instead of lower_bound_angle to deminish the number of coronas
    ttriangles_q = tripl(q)
    #min_angles_q = [(tight_angle(t)+3*lower_bound_angle(*t))/4 for t in ttriangles_q] # for all
    #if q==1: # for 26,23,35
    #    min_angles_q = [1.5*tight_angle(t) for t in ttriangles_q] # !!! for 26
    #else:
    min_angles_q = [(tight_angle(t)+lower_bound_angle(*t))/2 for t in ttriangles_q] # !!! for 144
    for vec in odo(min_angles_q, 2*RIFpi):
        if coronable(vec):
            res = []
            for (triangle, multiplicity) in zip(ttriangles_q, vec):
                res+=[triangle]*multiplicity
            yield (res)

def only_tight_corona_generator(q):
    # using a greater value instead of lower_bound_angle to deminish the number of coronas
    ttriangles_q = tripl(q)
    min_angles_q = [tight_angle(t) for t in ttriangles_q] # !!! for 26
    for vec in odo(min_angles_q, 2*RIFpi):
        if coronable(vec):
            res = []
            for (triangle, multiplicity) in zip(ttriangles_q, vec):
                res+=[triangle]*multiplicity
            yield (res)
            

#-------------------------------almost coronas
def find_almost_coronas(q):
    small_number=0.4
    almost_coronas=[]
    num = 0
    for triangles in corona_generator(q):
        sum_v = sum(V(*t) for t in triangles)
        sum_angles = sum(tight_angle(t) for t in triangles)
        diff = abs(2*RIFpi - sum_angles)
        if abs(diff)<small_number:
            if not diff.contains_zero():
                almost_coronas+=[[triangle_name(tt) for tt in triangles]]
                #draw_corona(q, corona_s)
            else:
                print("contains ZERO---------------------")
            print (n(diff))
            print(sum_v)
            print([triangle_name(tt) for tt in triangles])
            num+=1
    return(almost_coronas)
            
#------------------------------- For testing

def Q_inequality(x,y,z):
    r_min, r_max = min(x,y,z), (x+y+z-min(x,y,z))/2
    Q = LQ[(x,y)][1] + LQ[(x,z)][1] + LQ[(z,y)][1]
    return Q <= 1/s*(1+r_max/(r_max+s))*RIFpi*(r_min**2)/2


