#============================== Import from the table ==============================

ttriangles_symb = ['111', '11r', '11s', '1rr', '1rs', '1ss', 'rrr', 'rrs', 'rss', 'sss']

def read_table():
    '''Given a txt-table named name, returns its head and its rows as
    lists'''
    rows = []
    name="/home/tooticki/Research/DiscPackings/code/counter_examples_binary/input/cases.csv"
    with open(name) as f:
        s=(f.readline())[:-1] #remove \n
        head = (s).split(",")
        for line in f:
            if '0'<=line[0]<='9':
                line = line[:-1] #remove \n
                rows.append(line.split(","))
    return (head, rows)

def rows_to_dicts(head, rows):
    '''Given the head and the rows of a table, returns a list of
    dicts with adabted types of values'''
    dicts = []
    for row in rows:
        d = { h:r for (h,r) in zip(head, row)}
        for x in ['density', 'r', 's']:
            if d[x]!='':
                d[x] = float(d[x])  
        for x in ['case', '1-discs', 'r-discs', 's-discs', '111', '11r', '11s', '1rr', '1rs', '1ss', 'rrr', 'rrs', 'rss', 'sss']:
            if d[x]!='':
                d[x] = int(d[x])        
        dicts.append(d)
    return dicts

def row_to_dict(head, row):
    '''Given the head and the row  of a table, returns a
    dict with adabted types of values'''
    d = { h:r for (h,r) in zip(head, row)}
    for x in ['density', 'r', 's']:
        if d[x]!='':
            d[x] = float(d[x])  
    for x in ['case', '1-discs', 'r-discs', 's-discs', '111', '11r', '11s', '1rr', '1rs', '1ss', 'rrr', 'rrs', 'rss', 'sss']:
        if d[x]!='':
            d[x] = int(d[x])        
    return d


def input_values(case):
   (h,r) = read_table()
   d = row_to_dict(h,r[case-1])
   minpoly_r = (d["minpoly r"]).replace('r', 'x')
   minpoly_s = (d["minpoly s"]).replace('s', 'x')
   var('x')
   r_list = PolynomialRing(QQ,x)(minpoly_r).roots(AA,multiplicities=false)
   s_list = PolynomialRing(QQ,x)(minpoly_s).roots(AA,multiplicities=false)
   r = find_closest_root(r_list, d['r'])
   s = find_closest_root(s_list, d['s'])
   disc_nums_1rs =   [d[dt] for dt in ['1-discs', 'r-discs', 's-discs']]
   triangles_nums =  [d[t] for t in ttriangles_symb]
   coronas = [d['s-corona'],d['r-corona']]
   dens = idensity(*([r,s]+ disc_nums_1rs+triangles_nums))
   return(RIF(r),RIF(s), disc_nums_1rs, triangles_nums, coronas, dens)
  


def idensity(r, s, n1, nr, ns, m111, m11r, m11s, m1rr, m1rs, m1ss, mrrr, mrrs, mrss, msss):
    '''Density of the packing whose domain consists of nx circles of
    radius x and mxyz triangles formed by ciclres of radii x,y,z'''
    cov = pi*(n1 + nr*(r**2) + ns*(s**2))
    triangles = m111*[[1,1,1]]+m11r*[[1,1,r]]+m11s*[[1,1,s]]+m1rr*[[1,r,r]]+m1rs*[[1,r,s]]+m1ss*[[1,s,s]]+mrrr*[[r,r,r]]+mrrs*[[r,r,s]]+mrss*[[r,s,s]]+msss*[[s,s,s]]
    area = sum([d_area(*t) for t in triangles])
    return cov / area
