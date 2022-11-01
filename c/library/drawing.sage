
# ------------------------------------------------ Drawing ------------------------

def find_center(r, neighbors):
    '''Returns the coordinates of a new disc with two neighbors given as a
    pair of dictionaries  '''
    (x1,y1), (x2,y2) = neighbors[0]['xy'], neighbors[1]['xy']
    dx, dy = x2-x1, y2-y1
    r1, r2  =  neighbors[0]['r'], neighbors[1]['r']
    a1, ratio = angle(r+r2,r+r1,r1+r2), (r1+r)/(r1+r2)
    xv, yv = turn(ratio*dx, ratio*dy, a1)
    center = (x1+xv, y1+yv)
    return(center)

def draw_packing(p_type,d,discs,shifts,a_range):
    '''Gets a list of discts describing the packing and returns a figure
    representing the domain of this packing'''
    a_range = [n(x) for x in a_range]
    if p_type == 'bin':
        eqsymb = "\\approx "
        dens = n(d,digits = 6)
    elif p_type == 'ter':
        eqsymb = "\\leq "
        dens = n(d.upper(),digits = 6)
    else:
        eqsymb = "\\geq "
        dens = n(RIF(d).lower(),digits = 6)
    xmin,xmax,ymin,ymax = tuple(a_range)
    col = lambda x: (227/256,15/256,15/256) if x=='s' else (255/256,212/256,40/256) if x=='1' else (84/256,151/256,255/256)
    figure = Graphics()
    for shift in shifts:
        for disc in discs:
            figure += circle(vector(disc['xy'])+vector(shift), disc['r'], fill=True, facecolor = col(disc['rname']), edgecolor = 'black', thickness = .1, zorder=1)

    R,r = max([disc['r'] for disc in discs]) , min([disc['r'] for disc in discs])
    Rname,rname = min([disc['rname'] for disc in discs]), max([disc['rname'] for disc in discs])
    #rr = [x for x in [disc['r'] for disc in discs if disc['rname']=='r']][0] for ter!!!
    border_width = (xmax-xmin)/10
    white_rec = polygon([(xmin-border_width,ymin),(xmax+border_width,ymin),(xmax+border_width,ymin-border_width*2),(xmin-border_width,ymin-border_width*2)],color='white',zorder=1)
    figure+=white_rec
    figure += text(r'$\delta'+eqsymb+str(dens)+'$',((xmin+xmax)/2,ymin-border_width/1.5),horizontal_alignment = 'right',color='black',fontsize=23 )
    if Rname == '1':
        figure += text(r'$'+rname+'\\approx '+str(n(r,digits=6))+'$', ((3.5*xmin+4*xmax)/7.5,ymin-border_width/1.5),horizontal_alignment = 'left', color=col(rname),fontsize=23 )
    #figure += text("r="+str(n(rr,digits=8)), ((3*xmin+4*xmax)/7,ymin-6*border_width/7), horizontal_alignment = 'left', color=col('r'),fontsize='large' )
    else:
        figure += text(r'$\frac{'+rname+'}{'+Rname+"}\\approx "+str(n(r/R,digits=5))+'$', ((3.6*xmin+4*xmax)/7.6,ymin-border_width/1.5), horizontal_alignment = 'left', color=col(Rname),fontsize=23 )
    figure.set_axes_range(xmin,xmax,ymin-border_width,ymax)
    return figure



def baricenter(points):
    N = len(points)
    x = sum([p[0] for p in points])/N
    y = sum([p[1] for p in points])/N
    return(x,y)


def draw_quasi_b1_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    rr = R*(bin_dict['b1'][0])
    discs=[]
    discs+=[{'r':r, 'xy':(-r,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(r,0), 'rname':rname}]
    xc,yc = 0,0
    discs+=[{'r':R, 'xy':(0,sqrt((R+rr)**2-rr**2)), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(0,-sqrt((R+rr)**2-rr**2)), 'rname':Rname}]
    a = pi-arccos(rr/(rr+R))-2*arcsin(R/(R+rr))
    xa, ya = rr + (rr+R)*cos(a), (rr + R)*sin(a)
    discs+=[{'r':R, 'xy':(xa,ya), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(xa,-ya), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-xa,-ya), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-xa,ya), 'rname':Rname}]    
    return draw_packing(p_type,d,discs,[(0,0)],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b1_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':r, 'xy':(r,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-r,0), 'rname':rname}]
    xc,yc = 0,0
    x2,y2 = find_center(R,discs)
    discs+=[{'r':R, 'xy':(x2,y2), 'rname':Rname}]
    x3,y3 = find_center(R,[discs[0],discs[-1]])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    x4,y4 = find_center(R,[discs[0],discs[-1]])
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-x2,-y2), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-x3,-y3), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-x4,-y4), 'rname':Rname}]
    
    v1,v2 = vector((x4-x2,y4-y2)),vector((x2+x3,y2+y3))
    return draw_packing(p_type,d,discs,[([0,0]),v1,-v1,-v2,v2,v1+v2,-v1-v2,],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b3_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':R, 'xy':(0,0), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(0,2*R), 'rname':Rname}]
    x2,y2 = find_center(r,discs)
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    x3,y3 = find_center(r,[discs[1],discs[0]])
    discs+=[{'r':r, 'xy':(x3,y3), 'rname':rname}]
    xc,yc = ((x2+x3)/2,(y2+y3)/2)
    
    x4,y4 = (x2-R-r,y2)
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]

    x5,y5 =(x4,-y4)
    discs+=[{'r':R, 'xy':(x5,y5), 'rname':Rname}]
    
    x6,y6 = find_center(r,[discs[-2],discs[-1]])
    discs+=[{'r':r, 'xy':(x6,y6), 'rname':rname}]
    
    v1,v2 = vector((0,2*R)),vector((x6-x3,y6-y3))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b3_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':r, 'xy':(0,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(0,2*r), 'rname':rname}]
    xc,yc = (0,r)
    x2,y2 = find_center(R,discs)
    discs+=[{'r':R, 'xy':(x2,y2), 'rname':Rname}]
    x3,y3 = find_center(R,[discs[1],discs[0]])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    x4,y4 = find_center(r,[discs[2],discs[1]])
    discs+=[{'r':r, 'xy':(x4,y4), 'rname':rname}]

    v1,v2 = vector((x4,y4)),vector((x3-x2,y3-y2))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])


def draw_quasi_b4_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':R, 'xy':(0,0), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(2*R,0), 'rname':Rname}]

    x2,y2 = find_center(r,discs)
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    xc,yc = x2,y2

    x3,y3 = find_center(R,[discs[-1],discs[1]])
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]

    x4,y4 = x3-2*R,y3
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]

    v1,v2 = vector((2*R,0)),vector((x4,y4))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b4_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':R, 'xy':(-R,0), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(R,0), 'rname':Rname}]
    
    x2,y2 = find_center(r,discs)
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    xc,yc = x2,y2
    
    x3,y3 = find_center(R,[discs[0],discs[-1]])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    
    x4,y4 = find_center(R,[discs[-2],discs[1]])
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]

    x5,y5 = (0,2*y3-y2)
    discs+=[{'r':r, 'xy':(x5,y5), 'rname':rname}]

    discs+=[{'r':R, 'xy':(-R,2*y3), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(R,2*y3), 'rname':Rname}]

    v1,v2 = vector((x4+R,y4)),vector((x3-R,y3))
    
    return draw_packing(p_type,d,discs,[(0,0),v1,v2,-v1-v2,v1+v2,-v1,-v2],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])



def draw_quasi_b5_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':r, 'xy':(0,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(2*r,0), 'rname':rname}]
    xc,yc=0,0
    
    x2,y2 = find_center(r,discs)
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    x3,y3 = find_center(r,[discs[0],discs[-1]])
    discs+=[{'r':r, 'xy':(x3,y3), 'rname':rname}]
    x4,y4 = find_center(r,[discs[0],discs[-1]])
    discs+=[{'r':r, 'xy':(x4,y4), 'rname':rname}]
    x5,y5 = find_center(r,[discs[0],discs[-1]])
    discs+=[{'r':r, 'xy':(x5,y5), 'rname':rname}]
    x6,y6 = find_center(r,[discs[0],discs[-1]])
    discs+=[{'r':r, 'xy':(x6,y6), 'rname':rname}]

    x2,y2 = find_center(R,[discs[1],discs[6]])
    discs+=[{'r':R, 'xy':(x2,y2), 'rname':Rname}]
    x3,y3 = find_center(R,[discs[2],discs[1]])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    x4,y4 = find_center(R,[discs[3],discs[2]])
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]
    x5,y5 = find_center(R,[discs[4],discs[3]])
    discs+=[{'r':R, 'xy':(x5,y5), 'rname':Rname}]
    x6,y6 = find_center(R,[discs[5],discs[4]])
    discs+=[{'r':R, 'xy':(x6,y6), 'rname':Rname}]
    x7,y7 = find_center(R,[discs[6],discs[5]])
    discs+=[{'r':R, 'xy':(x7,y7), 'rname':Rname}]

    v1,v2 = vector((x5-x3,y5-y3)),vector((x6-x4,y6-y4))
    return draw_packing(p_type,d,discs,[(0,0),v1,v2,v1+v2,v1-v2,v2-v1,-v1-v2,-v1,-v2],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])


def draw_quasi_b6_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':r, 'xy':(-r,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(r,0), 'rname':rname}]
    x2,y2 = find_center(r,discs)
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    
    x3,y3 = find_center(R,[discs[1],discs[0]])
    xc,yc = x3,y3
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    x4,y4 = find_center(R,[discs[2],discs[1]])
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]
    x5,y5 = find_center(R,[discs[0],discs[2]])
    discs+=[{'r':R, 'xy':(x5,y5), 'rname':Rname}]
    
    x6,y6 = (0,2*y4-y2)
    discs+=[{'r':r, 'xy':(x6,y6), 'rname':rname}]
    x7,y7 = (-r,2*y4)
    discs+=[{'r':r, 'xy':(x7,y7), 'rname':rname}]
    x8,y8 = (r,2*y4)
    discs+=[{'r':r, 'xy':(x8,y8), 'rname':rname}]

    x9,y9 = find_center(R,[discs[7],discs[8]])
    discs+=[{'r':R, 'xy':(x9,y9), 'rname':Rname}]
    
    v1,v2 = vector((x5-x4,0)),vector((0,y9-y3))
    return draw_packing(p_type,d,discs,[(0,0)]+[(i*v1+j*v2+k*(v1+v2)/2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3] for k in [-1,0,1]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])



def draw_quasi_b6_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    XX = 2*sin(5*pi/12)*2*r
    xc,yc = (0,0)
    discs+=[{'r':r, 'xy':(0,XX), 'rname':rname}]
    discs+=[{'r':r, 'xy':(XX,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(0,-XX), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-XX,0), 'rname':rname}]

    x1,y1 = (sqrt(3)*XX/2,XX/2)
    
    discs+=[{'r':r, 'xy':(x1,y1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-x1,-y1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(y1,x1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-y1,-x1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(x1,-y1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-x1,y1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(y1,-x1), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-y1,x1), 'rname':rname}]

    discs+=[{'r':R, 'xy':(0,0), 'rname':Rname}]

    v1,v2 = vector((y1,-x1-XX)),vector((-x1-XX,y1))
    return draw_packing(p_type,d,discs,[(0,0)]+[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])


# TODO
def draw_quasi_b7_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    discs=[] 
    a = pi - tight_angle((R,r,R))
    discs+=[{'r':R, 'xy':(0,(r+R)*sin(a)), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(0,-(r+R)*sin(a)), 'rname':Rname}]
    xc,yc = 0,0
    discs+=[{'r':r, 'xy':((r+R)*cos(a),0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(-(r+R)*cos(a),0), 'rname':rname}]
    discs+=[{'r':R, 'xy':(-(r+R)*cos(a)-r-R,0), 'rname':Rname}]
    discs+=[{'r':R, 'xy':((r+R)*cos(a)+r+R,0), 'rname':Rname}]
    x1,y1 = find_center(r,[discs[0],discs[-1]])
    discs+=[{'r':r, 'xy':(x1,y1), 'rname':rname}]
    
    v1,v2 = vector((2*((r+R)*cos(a)+r+R), 0)),vector((x1+(r+R)*cos(a),y1))
    return draw_packing(p_type,d,discs,[(0,0)]+[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3] for k in [-1,0,1]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b7_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    discs=[]
    discs+=[{'r':r, 'xy':(-r,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(r,0), 'rname':rname}]
    xc,yc = 0,0
    x2,y2 = find_center(R,discs)
    discs+=[{'r':R, 'xy':(x2,y2), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(x2,-y2), 'rname':Rname}]
    x3,y3 =  find_center(R,[discs[0],discs[2]])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-x3,-y3), 'rname':Rname}]
    v1,v2 = vector((0,2*y2)),vector((2*x3,2*y3))                                   
    return draw_packing(p_type,d,discs,[(i*v1+j*v2+k*(v1+v2)/2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3] for k in [-1,0,1]], [xc-3*R,xc+3*R,yc-3*R,yc+3*R])

def draw_quasi_b8_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    disc0={'r':R, 'xy':(0,0), 'rname':Rname}
    disc1={'r':R, 'xy':(2*R,0), 'rname':Rname}
    discs = [disc0,disc1]
    discs+=[{'r':R, 'xy':find_center(R,[disc0,disc1]), 'rname':Rname}]
    x2,y2 = baricenter([dd['xy'] for dd in discs])
    discs+=[{'r':r, 'xy':(x2,y2), 'rname':rname}]
    discs+=[{'r':r, 'xy':(x2,-y2), 'rname':rname}]
    xc,yc = x2,y2
    v1,v2 = vector((2*R,0)),vector(find_center(R,[disc0,disc1]))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])
    
def draw_quasi_b8_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    disc0={'r':r, 'xy':(0,0), 'rname':rname}
    disc1={'r':R, 'xy':(0,r+R), 'rname':Rname}
    discs = [disc0,disc1]
    x2,y2 = find_center(R,discs)
    discs+=[{'r':R, 'xy':(x2,y2), 'rname':Rname}]
    discs+=[{'r':R, 'xy':(-x2,y2), 'rname':Rname}]
    discs+=[{'r':r, 'xy':(0,2*y2), 'rname':rname}]
    xc,yc = (0,0)
    discs+=[{'r':R, 'xy':(0,2*y2-r-R), 'rname':Rname}]
    v1,v2 = vector((2*x2,0)),vector((0,2*y2-2*r-2*R))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2+k*(v1+v2)/2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3] for k in [-1,0,1]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])


def draw_quasi_b9_r_smaller(d,r,R,rnames):
    p_type = "cex"
    rname,Rname = rnames[0],rnames[1]
    disc0={'r':R, 'xy':(0,0), 'rname':Rname}
    disc1={'r':R, 'xy':(2*R,0), 'rname':Rname}
    discs = [disc0,disc1]
    discs+=[{'r':R, 'xy':find_center(R,[disc0,disc1]), 'rname':Rname}]
    xc,yc = baricenter([di['xy'] for di in discs])
    discs+=[{'r':r, 'xy':(xc+r,yc+r/sqrt(3)), 'rname':rname}]
    discs+=[{'r':r, 'xy':(xc-r,yc+r/sqrt(3)), 'rname':rname}]
    discs+=[{'r':r, 'xy':(xc,yc-2*r/sqrt(3)), 'rname':rname}]
    discs+=[{'r':r, 'xy':(xc+r,-yc-r/sqrt(3)), 'rname':rname}]
    discs+=[{'r':r, 'xy':(xc-r,-yc-r/sqrt(3)), 'rname':rname}]
    discs+=[{'r':r, 'xy':(xc,-yc+2*r/sqrt(3)), 'rname':rname}]
    v1,v2 = vector((2*R,0)),vector(find_center(R,[disc0,disc1]))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])
    
def draw_quasi_b9_r_greater(p_type,d,r,R,rnames):
    rname,Rname = rnames[0],rnames[1]
    disc0={'r':r, 'xy':(0,0), 'rname':rname}
    disc1={'r':r, 'xy':(0,2*r), 'rname':rname}
    discs = [disc0,disc1]
    x2,y2 = find_center(r,discs)
    disc2={'r':r, 'xy':(x2,y2), 'rname':rname}
    xc,yc = baricenter([di['xy'] for di in discs])
    discs+=[disc2]
    x3,y3 = find_center(R,[disc1,disc0])
    discs+=[{'r':R, 'xy':(x3,y3), 'rname':Rname}]
    x4,y4 = find_center(R,[disc2,disc1])
    discs+=[{'r':R, 'xy':(x4,y4), 'rname':Rname}]
    x5,y5 = find_center(R,[disc0,disc2])
    discs+=[{'r':R, 'xy':(x5,y5), 'rname':Rname}]
    
    discs+=[{'r':r, 'xy':(2*x4,0), 'rname':rname}]
    discs+=[{'r':r, 'xy':(2*x4,2*r), 'rname':rname}]
    discs+=[{'r':r, 'xy':find_center(r,[discs[-1],discs[-2]]), 'rname':rname}]
    
    v1,v2 = vector((x5-x3,y5-y3)),vector((x4-x3,y4-y3))
    return draw_packing(p_type,d,discs,[(i*v1+j*v2) for i in [-3,-2,-1,0,1,2,3] for j in [-3,-2,-1,0,1,2,3]],[xc-3*R,xc+3*R,yc-3*R,yc+3*R])


def draw_and_save_found_cases():
    for c in found_saturated:
        cd = best_bin_for_found[c]
        binname,rnames = cd['bin_name'],cd['rtype']
        d_new,dter= cd['d_new'],cd['ter_d']
        ineq,r,R = cd['sign'],cd['r'],cd['R']
        dter = ter_dict[c][-1]
        f = Graphics()
        if binname =='b1':
            if ineq == '<':
                f = draw_quasi_b1_r_smaller(d_new,r,R,rnames)               
            else:
                f = draw_quasi_b1_r_greater("cex",d_new,r,R,rnames)               
        if binname =='b3':
            if ineq == '<':
                f = draw_quasi_b3_r_smaller(d_new,r,R,rnames)               
            else:
                f = draw_quasi_b3_r_greater("cex",d_new,r,R,rnames)               
        if binname =='b4':
            if ineq == '<':
                f = draw_quasi_b4_r_smaller(d_new,r,R,rnames)               
            else:
                f = draw_quasi_b4_r_greater("cex",d_new,r,R,rnames)
        if binname =='b5':
            if ineq == '>':
                f = draw_quasi_b5_r_greater("cex",d_new,r,R,rnames)
        if binname =='b6':
            if ineq == '>':
                f = draw_quasi_b6_r_greater("cex",d_new,r,R,rnames)
            else:
                f = draw_quasi_b6_r_smaller(d_new,r,R,rnames)
        elif binname == 'b7':
            if ineq == '<':
                f = draw_quasi_b7_r_smaller(d_new,r,R,rnames)               
            else:
                f = draw_quasi_b7_r_greater("cex",d_new,r,R,rnames)               
        elif binname == 'b8':
            if ineq == '<':
                f = draw_quasi_b8_r_smaller(d_new,r,R,rnames)                
            else:
                f = draw_quasi_b8_r_greater("cex",d_new,r,R,rnames)
        elif binname == 'b9':
            if ineq == '<':
                f = draw_quasi_b9_r_smaller(d_new,r,R,rnames)
            else:
                f = draw_quasi_b9_r_greater("cex",d_new,r,R,rnames)
        #f.save(dpi=800, axes = False, filename = "pictures/"+str(c)+"_"+binname+"_r"+ineq+"_"+rnames+".pdf")
        f.save(dpi=100, axes = False, filename = "pictures/"+str(c)+"_counter.pdf")





#============================== Draw Triangulated Ternary Packings ============================================================


def encoding(case,r,s):
    '''Returns the encoding of the case case where the radii are equal to
    r and s'''
    from packings3discs import packings
    return (packings(r,s))[case-1]

def halfplane_sign(p, a, b):
    '''Returns the sign of the point p with respect to the line passing
    by a and b'''
    px,py = p
    ax,ay = a
    bx,by = b
    # y = mx + n  is the line of the points a and b
    m = (by-ay)/(bx-ax)
    n = ay-m*ax
    return (py - (m*px + n)).sign()


def find_center(r, neighbors):
    '''Returns the coordinates of a new disc with given neighbors as a
    pair of dictionaries'''
    (x1,y1), (x2,y2) = neighbors[0]['xy'], neighbors[1]['xy']
    dx, dy = x2-x1, y2-y1
    r1, r2  =  neighbors[0]['r'], neighbors[1]['r']
    a1, ratio = angle(r+r2,r+r1,r1+r2), (r1+r)/(r1+r2)
    xv, yv = turn(ratio*dx, ratio*dy, a1)
    center = (x1+xv, y1+yv)
    return(center)

def touch(r1, x1, y1, r2, x2, y2, epsilon):
    return sqrt((x1-x2)^2 + (y1-y2)^2) - (r1+r2) < epsilon

# TODO: test it better
def other_neighbors(r, xy, neighbors, discs, epsilon):
    '''Returns the neighbors of a new disc except the already given ones
    with given precision epsilon '''
    on = []
    n_discs = [d for d in discs if not d in neighbors]
    for disc in n_discs :
        if touch(r, xy[0], xy[1], disc['r'], disc['xy'][0], disc['xy'][1], epsilon):
            on.append(discs.index(disc))
    #if len(on) > 0:
    #    print("Other neighbors: "+str(len(on)))
    return(on)

# TODO: test better ?
def find_epsilon(r,s):
    '''Returns the mimimal distanc between non-tangent discs in the
    compact packing'''
    radii = [1,r,s]
    pairs = [(r1,r2) for r1 in radii for r2 in radii if r1<=r2]
    min_dist = 2
    for (r_top,r_bottom) in pairs:
        for (r_left, r_right) in pairs:
            a_lt = angle(r_left+r_bottom, r_left+r_top, r_top+r_bottom)
            a_rt = angle(r_left+r_bottom, r_left+r_top, r_top+r_bottom)
            dist = sin(a_lt)*(r_left+r_top) + sin(a_rt)*(r_right+r_top) - r_left -r_right
            if dist < min_dist:
                min_dist = dist
    return min_dist

# turn (x,y) by 'a' counterclockwise
turn = lambda x,y,a: (x*cos(a)-y*sin(a), x*sin(a)+y*cos(a))

# OK
def discs(case,r,s):    
    '''Returns the geometric properties of the packing'''
    l = (encoding(case,r,s))[1]
    name = lambda x : '1' if x>r else 'r' if x>s else 's'
    disc_0 = {'r':RIF(l[0]), 'xy':(0,0),'rname':name(l[0])}
    disc_1 = {'r':RIF(l[1]), 'xy':(disc_0['r']+RIF(l[1]),0), 'rname':name(l[1])}
    discs = [disc_0,disc_1]
    # the precision: two discs can not be at distance < epsilon and do not touch
    epsilon = find_epsilon(r,s)/2
    for entry in l[2:]:
        R = RIF(entry[2])
        neighbors = list(entry[:2])
        xy = find_center(R, [discs[i] for i in neighbors])
        d = {'r':R, 'xy':xy, 'rname':name(R)}
        index = l.index(entry)
        if len(entry)==4:
            twin = entry[3]
            (discs[twin]).update({'twin':index})
            d.update({'twin':twin})
        discs.append(d)
        #print("Disc number "+str(index))# + ":   " + str([r,xy,neighbors[:2]]))


    periods = []
    for disc in discs:
        if 'twin' in disc.keys():
            disctwin = discs[disc['twin']]
            x1,y1 = disc['xy']
            x2,y2 = disctwin['xy']
            v = vector((x1-x2,y1-y2))
            if (not v in periods) and (not (-v) in periods) :
                periods.append(v)
    return discs,periods



def save_found_ter():
    for case in found_saturated:
        r,s,sr,d = ter_dict[case]
        all_discs,periods = discs(case,r,s)
        binname = best_bin_for_found[case]['bin_name']
        shifts = [i*v1+j*v2 for i in [-2,-1,0,1,2] for j in [-2,-1,0,1,2] for v1 in periods for v2 in periods]
        f = draw_packing("ter",d,all_discs,shifts,[-3,3,-3,3])
        f.save(dpi=100, axes = False, filename = "pictures/"+str(case)+"_ter.pdf")

def save_one_ter(case):    
        r,s,sr,d = ter_dict[case]
        all_discs,periods = discs(case,r,s)
        binname = best_bin_for_found[case]['bin_name']
        shifts = [i*v1+j*v2 for i in [-2,-1,0,1,2] for j in [-2,-1,0,1,2] for v1 in periods for v2 in periods]
        f = draw_packing("ter",d,all_discs,shifts,[-3,3,-3,3])
        f.save(dpi=100, axes = False, filename = "pictures/"+str(case)+"_one_ter.pdf")



def save_all_bin():
    for binname in bin_names:
        d,r,R,rnames = bin_dict[binname][1],bin_dict[binname][0],1,['r','1']
        f = Graphics()
        if binname =='b1':
            f = draw_quasi_b1_r_greater("bin",d,r,R,rnames)               
        elif binname =='b3':
            f = draw_quasi_b3_r_greater("bin",d,r,R,rnames)               
        elif binname =='b4':
            f = draw_quasi_b4_r_greater("bin",d,r,R,rnames)
        elif binname =='b5':
            f = draw_quasi_b5_r_greater("bin",d,r,R,rnames)
        elif binname =='b6':
            f = draw_quasi_b6_r_greater("bin",d,r,R,rnames)
        elif binname == 'b7':
            f = draw_quasi_b7_r_greater("bin",d,r,R,rnames)               
        elif binname == 'b8':
            f = draw_quasi_b8_r_greater("bin",d,r,R,rnames)
        elif binname == 'b9':
            f = draw_quasi_b9_r_greater("bin",d,r,R,rnames)
        f.save(dpi=100, axes = False, filename = "pictures/"+binname+".pdf")
        
