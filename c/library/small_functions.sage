#------------------------------ Small Useful Functions ------------------------------

def stick(xyr1,xyr2,r3,direct=true):
    (x1,y1,r1),(x2,y2,r2) =  xyr1,xyr2
    l=sqrt((x2-x1)^2+(y2-y1)^2) # l=r1+r2 if discs 1 & 2 are tangents (not mandatory)
    cosa=((r1+r3)^2+l^2-(r2+r3)^2)/(2*(r1+r3)*l)
    sina=(1 if direct else -1)*sqrt(1-cosa^2)
    x3=x1+(cosa*(x2-x1)-sina*(y2-y1))*(r1+r3)/l
    y3=y1+(sina*(x2-x1)+cosa*(y2-y1))*(r1+r3)/l
    return (x3.canonicalize_radical(),y3.canonicalize_radical(),r3)


# find the algebraic roots in [a,b] of a radical expression f
def find_roots(f,a,b):
    f=numerator(f.canonicalize_radical())
    g=expand(f)
    r=g.variables()[0]
    while not g.is_polynomial(r):
        y=choose_sqrt_in(g) # choose a sqrt to be removed
        x=g.coefficient(sqrt(y))
        z=g-x*sqrt(y) # g=x*sqrt(y)+z
        g=x^2*y-z^2
        g=g.canonicalize_radical()
    g=AA[r](g) # QQ
    return [i for i in g.roots(AA,multiplicities=false) if i>=a and i<=b and RIF(f.subs(r=i)).contains_zero()]



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

# radius value by name
r_value = lambda x:1 if x=='1' else r if x=='r' else s

#radius name by value
r_name = lambda x:'1' if x>r else 'r' if x>s else 's'

# triangle type detection
triangle_name = lambda t: r_name(t[0])+r_name(t[1])+r_name(t[2])

# area of a triangle with sides of length a, b and c
area = lambda a,b,c: sqrt((a+b+c)*(a+b-c)*(a-b+c)*(b+c-a))/4

# area of a tight triangle with discs of radii ra, rb and rc
d_area = lambda ra,rb,rc: area(ra+rb,rb+rc,ra+rc)

# angle opposite to edge a
angle = lambda a,b,c:arccos((b^2+c^2-a^2)/(2*b*c))

# given t=(ra,rb,rc) returns the angle in the center of the circle rb
tight_angle = lambda t: RIF(angle(t[0]+t[2], t[0]+t[1], t[1]+t[2]))


# area of the triangle with sides of length a, b and c covered by the discs of radius ra, rb and rc
cov = lambda a,b,c,ra,rb,rc: angle(a,b,c)*ra^2/2+angle(b,c,a)*rb^2/2+angle(c,a,b)*rc^2/2


# given q, returns all triples of radii with q in the middle
tripl = lambda q: [(a,q,b) for (a,b) in [(1,1),(1,r),(1,s),(r,r),(r,s),(s,s)]]

