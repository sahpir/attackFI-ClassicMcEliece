#!/usr/bin/env sage

"""
Author: Julian Danner
Descriptions: Sabine Pircher
Date: July 2022
This program runs currently only on Linux, because it implements a function call to an external C-program.
"""
import time
import argparse
import itertools
import os
import tempfile

#get number of cores
NCORE = len(os.sched_getaffinity(0))

#collection of funcs for solving linear and quadratic fault eqs and computing an alternative secret key
"""
Reads polynomials from file with this format: "%u %u %u %u\n"
Each polynomial is stored in a new line
"""
def read_polys(fname:str) -> 'list[list[int]]':
  polys = []
  with open(fname) as f:
    for line in f:
      polys.append( list( int(i) for i in line.split(' ') ) )
  return polys

"""
Reads polynomials from file with this format: "%u %u,%x\n" e.g. 1 3487,000111001100
Each polynomial is stored in a new line
"""
def read_lin_polys_with_fault(fname:str) -> 'list[tuple[list[int],str]]':
  polys = []
  with open(fname) as f:
    for line in f:
      line_ = line.strip().split(',')
      if len(line_)==1:
        polys.append( (list( int(i) for i in line_[0].split(' ') ),'') )
      elif len(line_)==2:
        polys.append( (list( int(i) for i in line_[0].split(' ') ),line_[1]) )
  return polys

"""
Reads codewords from file in binary format or hex format
"""
def read_code_words(fname:str) -> 'list[str]':
  #decide if repr is binary or hex! (hex must have spaces between hex-vals!)
  with open(fname) as f:
    first_line = f.readline()
  if first_line[2]==' ':
    return read_code_words_hex(fname)
  else:
    return read_code_words_bin(fname)

"""
Reads codewords from file in the format of 1 bytes, each byte is seperated by a space e.g. df f6 19 ...
Each codeword is stored in a new line.
"""
def read_code_words_hex(fname:str) -> 'list[str]':
  code_words = []
  with open(fname) as f:
    for line in f:
      line_bin = ''
      for h in line.strip().split(' '):
        h_bin = format(int(h, base=16), "08b")
        line_bin += h_bin[::-1]
      code_words.append( line_bin )
  return code_words

"""
Reads codewords from file in binary format eg. 01110011...
Each codeword is stored in a new line.
"""
def read_code_words_bin(fname:str) -> 'list[str]':
  cws = []
  with open(fname) as f:
    for line in f:
      cws.append( line.strip() )
  return cws

"""
Calculates the corresponding generator matrix from given public key in binary format
"""
def code_gen_from_pk(fname_pk:str, n:int, m:int, t:int):
  pk = []
  with open(fname_pk) as f:
    for line in f:
      pk.append( line.strip() )
  #assert( len(pk)==1 )
  pk = pk[0]
  #parse into matrix
  F2 = GF(2)
  H_ = Matrix(F2, m*t, n-m*t, [F2(e) for e in pk], sparse = False)
  #H = block_matrix( [1, H_], nrows=1 )
  #get generator matrix
  G = block_matrix( [H_.transpose(), 1], nrows=1 )
  #echelonize
  G.echelonize()
  return G

"""
Calculates the corresponding generator matrix from given public key that is formatted in hex with each byte sperated by a ","
"""
def code_gen_from_pk_hex(fname_pk_hex:str, n:int, m:int, t:int):
  with open(fname_pk_hex) as f:
    pk = ''.join( format(int(h, base=16), "08b")[::-1] for h in f.readline().split(',') if len(h)>0)
  #parse into matrix
  F2 = GF(2)
  H_ = Matrix(F2, m*t, n-m*t, [F2(e) for e in pk], sparse = False)
  #H = block_matrix( [1, H_], nrows=1 )
  #get generator matrix
  G = block_matrix( [H_.transpose(), 1], nrows=1 )
  #echelonize
  G.echelonize()
  return G


"""
Reads the support elements from file formatted in bytes sepearated by comma ',' 
Input: filename fname to read from, extension degree m
Return: list of support elements in F_{2^m}
"""
def read_alpha(m:int, fname:str) -> 'list[list[int]]':
  F2m=GF(2^m, names='t', modulus='minimal_weight') #todo is it constructed w.r.t. the correct modpoly?
  code_words = []
  with open(fname) as f:
    support_hex = f.readlines()[0].strip().split(',')
    #support_hex = f.readlines()[0].strip().split(' ')
  t = F2m('t')
  alpha_bin = [ bin(int(h, base=16))[2:] for h in support_hex ]
  alpha = [ sum([t**i *F2m(a[-(i+1)]) for i in range(len(a))], F2m.zero()) for a in alpha_bin ]
  return alpha

"""
Reads the Goppa Polynomial from file. Coefficients are formatted in bytes separated by comma ','
Input: filename fname to read from, extension degree m
Return: Goppa polynomial in F_{2^m}[x]
"""
def read_goppa_poly(m:int, fname:str):
  coeffs = read_alpha(m, fname)
  F2m = coeffs[0].parent()
  F2mx.<x> = PolynomialRing(F2m, 'x')
  g = sum( x**i * F2mx(coeffs[i]) for i in range(len(coeffs)) )
  return g

"""
Converts a element in F_{2^m} to hex-Format
"""
def to_hex(el):
  return "%0.4x" % el.integer_representation()

"""
Writes the full support set to given file
"""
def write_alpha(a:list, fname:str):
  with open(fname, 'w') as f:
    f.write( ','.join( to_hex(ai) for ai in a ) )
  return

"""
Writes the Goppa polynomial to given file
"""
def write_poly(g, fname:str):
  write_alpha(g.list(), fname)

"""
Computes the generator matrix from given support set and Goppa polynomial
Input: full support set a, Goppa polynomial g
Return: generator matrix G
"""
def comp_goppa_code_gen(a:list, g):
  C = codes.GoppaCode(g,a)
  G = C.generator_matrix()
  return G

"""
Computes codewords
uses Zassenhaus alg to compute intersection of vector spaces
If all indeterminates of the unknown support elements are contained in the polynomial equations, 
then zero_idxs and one_idxs are empty and codewords are directly read from generator matrix
Input: generator matrix G, zero indices zero_idxs, one indices one_idxs
Return: codewords cws of the Goppa code
"""
def get_code_words(G, zero_idxs:'list[int]', one_idxs:'list[int]'):
  #print(f'  constructing codewords with zero_idxs={str(zero_idxs)} and one_idxs={str(one_idxs)}')
  n = G.ncols()
  if len(zero_idxs)==0 and len(one_idxs)==0:
    return [ G.row(i) for i in range(G.nrows()) ]
  #constuct vector space containing all vecs with 0s at zero_idxs and 1s at one_idxs
  free_idxs = list( set(range(n)) - set(list(zero_idxs)+list(one_idxs)) )
  A = Matrix(GF(2), len(free_idxs), n)
  #assert( A==0 )

  #fill A
  for i in range(A.nrows()):
    k = free_idxs.pop()
    A[i,k] = 1
    for j in one_idxs:
      A[i,j] = 1
  #assert( len(free_idxs)==0 )

  #create blockmat 
  B = block_matrix( [[G,G],[A,0]] )
  B.echelonize()
  #find rank of left half!
  q = next( i for i in range(G.rank(),B.nrows()) if B[i,:n]==0 )
  q = next( i for i in range(G.nrows(),B.nrows()) if B[i,:n]==0 )
  #extract basis of intersection of row-space of G and row-space of A:
  G_int_A = B[q:,n:]
  cws = list( G_int_A.row(i) for i in range( G_int_A.rank() ) )
  #assert( all( all( c[i].is_zero() for i in zero_idxs ) for c in cws ) )
  if len(one_idxs)==0:
    return cws
  cw_basis_idxs = [ i for i in range(len(cws)) if all( cws[i][j].is_one() for j in one_idxs ) ]
  #assert( len(cw_basis_idxs)>0 )
  cws = [ cws[i]+cws[cw_basis_idxs[0]] for i in range(len(cws)) if not i in cw_basis_idxs ]
  #assert( all( all( c[i].is_one() for i in one_idxs ) for c in cws ) )
  #assert( len(cws)>0 )
  return cws

"""
checks whether linear polys hold for alpha
"""
def check_linear(alpha, lin_polys) -> bool:
  for l in lin_polys:
    if not sum( alpha[i] for i in l ).is_zero():
      return false
  return true


def check_quad(alpha, quad_polys) -> bool:
  for q in quad_polys:
    if not (alpha[q[0]]*alpha[q[1]] + alpha[q[2]]*alpha[q[3]]).is_zero():
      return false
  return true

"""
Solves linear equations 
Sets up a matrix with each row corresponding to one linear equations and calculates the echelon form of this matrix.
The remaining non-zero rows (<=n)  are returned.
Input: linear polynomial equations lin_polys, length n
Return: reduced linear polynomials
"""
def solve_linear(n:int, lin_polys:list) -> 'list[list[int]]':
  start_time = time.time()
  print('start reducing linear polys.')
  s_time = time.time()
  F2 = GF(2)
  ##create dense matrix within M4RI -- very fast reduction!
  A = matrix(F2, len(lin_polys), n, sparse=False);
  for i in range(len(lin_polys)):
    l = set(lin_polys[i])
    if len(l)==1 or len(l)==4 or len(l)==3: #make sure to not parse [1 2 1 2] as a1+a2 (!)
      for j in l:
        A[i,j] = F2.one()
  print(f'  created matrix (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  A.echelonize()
  print(f'  computed rref (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  A = A[:A.rank(),:]
  #read out simplified system:
  lin_polys_red = []
  strA = str(A).replace(' ','').split('\n')
  for r in strA:
    l = [i for i,e in enumerate(r[1:-1]) if e=='1']
    if len(l)>0:
      lin_polys_red.append( l )
  print(f'  finalized output (took {time.time()-s_time:.2f}s)')
  print(f'end reducing linear polys (took {time.time()-start_time:.2f}s)')
  return lin_polys_red

"""
Solves linear equations with the fault value explicitly given, faults generated by BIT SET (bitwise OR)
"""
#note injected fault is NOT xored but ORed to values, i.e., if ai = aj OR eps,
#     then we have:    ai[k] = (aj[k]+1)*(eps[k]+1) +1
#     or equivalently: eps[k] = 0 => ai[k] = aj[k]
#                      eps[k] = 1 => ai[k] = 1
#     this observation is used to construct a linear system in the 'coordinates'
#     of all alphas in order to find a (hopefully) unique solution!
#note if only the 0s or 1s of eps can be trusted, i.e., are reliable and should be used;
#     this function can be adapted very easily!
def solve_linear_wf(n:int, m:int, lin_polys:list) -> '':
  start_time = time.time()
  print('start reducing linear polys.')
  s_time = time.time()
  F2 = GF(2)
  ##create dense matrix within M4RI -- very fast reduction!
  A = matrix(F2, len(lin_polys)*m, n*m+1, sparse=False);
  for idx in range(len(lin_polys)):
    l = lin_polys[idx][0]
    eps = lin_polys[idx][1]
    if len(l)==1:
      i = l[0]
      for k in range(m):
        A[ m*idx+k, m*i+k] = F2.one()
    elif len(l)==2:
      i = l[0]
      j = l[1]
      assert(i!=j)
      for k in range(m):
        if eps[k]=='1':
          A[m*idx+k, m*j+k] = F2.one()
          A[m*idx+k, n*m] = F2.one()
          #we have a[j][k] = 1
        elif eps[k]=='0':
          A[m*idx+k, m*i+k] = F2.one()
          A[m*idx+k, m*j+k] = F2.one()
          #we have a[i][k] = a[j][k]
    else:
      raise ValueError('lin_polys with fault involve unexpected number of indices!')

  print(f'  created matrix (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  A.echelonize()
  print(f'  computed rref (took {time.time()-s_time:.2f}s)')

  assert( n*m not in A.pivots() )

  s_time = time.time()
  one  = []
  zero = []
  rem  = []
  #read out simplified system:
  for i in range( A.rank() ):
    if i%100==0: print(f'processing row {i} of {A.rank()}',end='\r')
    d = sorted(list(A.row(i).dict().keys()))
    if len(d)==1:
      zero.append( Integer(d[0]).quo_rem(m) )
    elif len(d)==2 and d[1] == n*m:
      one.append(  Integer(d[0]).quo_rem(m) )
    else:
      rem.append( [Integer(d_).quo_rem(m) for d_ in d] )
      #not completely determined!
  print(f'  finalized output (took {time.time()-s_time:.2f}s)')

  ##### debugging! #####
  det = set( one + zero )
  #check which alphas are determined completely!
  determined_as = [ a for a in range(n) if all( (a,k) in det for k in range(m) ) ]
  len(determined_as)

  #convert all alphas to string repr -- put 'x' when val is not determined
  alphas_str = [ ['x']*m for i in range(n) ]
  for i,k in zero:
    alphas_str[i][k] = '0'
  for i,k in one:
    alphas_str[i][k] = '1'
  alphas_str = [ ''.join(a) for a in alphas_str ]
  alphas_str_ = [] #stores tuples: str-repr of alpha, list of idxs (of alpha) that have this repr!
  idxs = set(range(n))
  while len(idxs)>0:
    i = idxs.pop()
    ai_str = alphas_str[i]
    I = [ j for j in idxs if alphas_str[i]==alphas_str[j] ]
    alphas_str_.append( (alphas_str[i], [i]+I) )
    idxs -= set(I)
  alphas_str_ = sorted(alphas_str_)
  tmp = [I for a,I in alphas_str_ if len(I)==a.count('x')]
  tmp = [(a,I) for a,I in alphas_str_ if len(I)==2**a.count('x')]
  #####

  print(f'end reducing linear polys (took {time.time()-start_time:.2f}s)')
  return lin_polys_red

"""
Solves the polynomial system of equations using the reduced linear equations and quadratic equations
Note: lin_polys must be reduced!
Determines the unknown support elements by solving the polynomial system of equations.
If the polynomial system of equations contain only a subset of indeterminates in the unknown support set,
 it outputs the determined subset of support elements and the remaining indices, the guessed indices and the updated indices.
If the polynomial system of equations contains all indeterminates in the unknown support set,
 rem_idxs, upd_idxs and guessed_idxs are empty
Input: reduced linear polynomials, quadratic polynomials
Return: support set S_L, indices of remaining unknowns rem_idxs, updated indices upd_idxs, guessed indices guessed_idx
"""
def solve_quad(n:int, m:int, linear_polys:list, quadratic_polys:list):
  start_time = time.time()
  print('start solving quadratic eqs.')
  s_time = time.time()
  F2m=GF(2^m, names='t', modulus='minimal_weight') #todo is it constructed w.r.t. the correct modpoly?
  ##comp number of remaining indets (note 1-indexed!)
  rem_inds = list( set(i for i in range(n)) - set(l[0] for l in linear_polys) )
  ##fix first rem_inds to 1, i.e., remove it
  one_idx = rem_inds[0]
  rem_inds = rem_inds[1:]
  ##create corr ring
  F2ma = PolynomialRing(F2m, names=['a'+str(k) for k in rem_inds])
  toF2m = lambda a : a.constant_coefficient() if sage.rings.polynomial.multi_polynomial_element.is_MPolynomial(a) and a.degree()==0 else a
  ##helper map from rem_inds to corr vars
  var_map = dict( zip(rem_inds,[F2ma('a'+str(k)) for k in rem_inds]) )
  ##fix one of remaining inds to 1
  var_map[ one_idx ] = F2ma.one()
  ##get linear reduction map
  l_red = {}
  l_red[one_idx] = F2m.one()
  for l in linear_polys:
    l_red[l[0]] = sum( [ var_map[l] for l in l[1:] ], F2m.zero() )
  ##complete l_red s.t. it maps all inds to their reductions
  for l in rem_inds:
    l_red[l] = var_map[l]
  ##create - reduced - quadratic eqs
  quad_eqs = []
  for q_idxs in quadratic_polys:
    q = l_red[ q_idxs[0] ] * l_red[ q_idxs[1] ] + l_red[ q_idxs[2] ]*l_red[ q_idxs[3] ]
    quad_eqs.append(q)
  print(f'  created reduced quadratic eqs (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  ##create ideal from quad_eqs
  I = F2ma.ideal(quad_eqs)
  I = F2ma.ideal( I.interreduced_basis() )
  print(f'  we have {len(rem_inds)} indets in {len(I.interreduced_basis())} equations')

  ##compute support candidate set
  #compute I_red -- only containing indets occurring in eqs
  F2ma_red = PolynomialRing(F2m, I.basis.variables())
  if F2ma.ngens()-F2ma_red.ngens()>0:
    print(f'  {F2ma.ngens()-F2ma_red.ngens()} of the remaining {F2ma.ngens()} indets did not occur in the quadratic polys; can only find partial support candidates!')
  I_red = I.change_ring(F2ma_red)
  d_red = I_red.dimension()
  #determine indices remaining - in case we can solve I_red (!)
  rem_idxs = set(rem_inds) - set( int(ind_name[1:]) for ind_name in F2ma_red.variable_names() )
  #determine incides that need to be updated after vals for rem_idxs are determined -- i.e., all idxs in which an index from rem_idxs occurs
  upd_idxs = set([i for i in range(n) if l_red[i].parent()!=F2m and any([ var_map[j] in l_red[i].variables() for j in rem_idxs ])]) - rem_idxs
  if d_red<=1:
    ind = I_red.ring().gens()[0]
    print(f'  using (truncated) exhaustive search on indet {ind} to shrink support candidate set...')
    #index that should be used for further processing -- since it was guessed!
    guessed_idxs = set( [int( str(ind)[1:] )] )
    #go through all F2m els (without 0 and 1 since they are already taken!)
    #TODO might lead to false results if index of 0 has not been found!
    L = set( F2m.list() ) - set( [F2m.zero(),F2m.one()] )
    #construct search space for ind, i.e., find G s.t. for every a in F2m/{0,1} there is g in G and j in {0,..,m-1} s.t. a = g^(2^j)
    G_ = []
    while len(L)>0:
      a = L.pop()
      G_.append( a )
      L -= set( a**(2**i) for i in range(m) )
    S_red_tmp = ( (I_red+F2ma_red.ideal(ind-a)).variety() for a in G_ )
    S_red = itertools.chain.from_iterable( S_red_tmp )
    print(f'  constructed generator of partial solutions (took {time.time()-s_time:.2f}s)           ')
    s_time = time.time()
    #extend to partial support candidates!
    S_red = ( { F2ma(k):s[k] for k in s } for s in S_red )
    S_L_ = ( list( toF2m(l_red[i].subs(sol)) for i in range(n) ) for sol in S_red )
    S_L = ( sol for sol in S_L_ if len(set(sol)) == len(sol) )
    if len(rem_idxs)>0:
      print(f'  constructed partial support candidate set (took {time.time()-s_time:.2f}s)')
    else:
      print(f'  constructed support candidate set (took {time.time()-s_time:.2f}s)')
  else:
    print(f'  reduced ideal has dimension {d_red} (>1); aborting computation -- restart with additional polynomials!')
    raise ValueError('support candidate set too large; i.e., too few equations!')
  
  print(f'end solving quadratic eqs (took {time.time()-start_time:.2f}s)')
  return S_L, rem_idxs, upd_idxs, guessed_idxs

"""
Computes the syndrome polynomial from support elements and non-zero positions of codeword 
Input: indetemindate x of F_{2^m}[x], support elements tildea, non-zero positions of codeword supp(c)
Return: syndrome polynomial in F_{2^m}[x]
"""
def comp_syn_poly(x, tildea:list, Ic:list):
  prod_tmp = product(x-tildea[i] for i in Ic)
  return sum( prod_tmp.quo_rem(x-tildea[i])[0] for i in Ic )

"""
Computes the square of the Goppa polynomial of degree 2*t from found support elements and given codewords. 
Codewords are generated by the generator matrix (not in this function).
Uses the Greatest Common Divider algorithm.
Input: generator matrix G, support set tildea, degree of Goppa polynomial t, codewords cws
Return: Goppa polynomial tildeg on success, 0 on failure, error if there are too less codewords to calculate tildeg 
"""
#note returns non-zero tildeg only if it was successful!
def goppaGCD(G:list, tildea:list, t:int, cws):
  #compute alternative goppa-polynomial - if there is one!
  start_time = time.time()
  print('start goppaGCD.')
  s_time = time.time()

  #construct ring
  F2m = tildea[0].parent()
  F2mx.<x> = PolynomialRing(F2m, 'x')
  tildeg = F2mx.zero()

  for c in cws:
    #syndrom poly can be computed!
    if tildeg.is_zero() or tildeg.degree()>2*t:
      Ic = list( i for i in range(len(c)) if c[i].is_one() )
      syn_poly = comp_syn_poly(x, tildea, Ic)
      tildeg = tildeg.gcd(syn_poly)
      print(f'  computed syndrome polynomial and gcd (took {time.time()-s_time:.2f}s), current degree {tildeg.degree()}')
      s_time = time.time()
    elif tildeg.degree()<2*t:
      print(f'end goppaGCD FAILURE, degree too low! (took {time.time()-start_time:.2f}s)')
      return F2mx.zero() #indicate failure!
    else:
      #degree is equal to 2*t
      print(f'end goppaGCD potential poly found! (took {time.time()-start_time:.2f}s)')
      return tildeg
  print(f'end goppaGCD FAILURE too few code-words supplied! (took {time.time()-start_time:.2f}s)')
  raise ValueError('start again with more code-words')

###BUGGY!
#def lin_inv_mod(lin, tildeg):
#  QR = tildeg.quo_rem( lin )
#  return QR[0]*QR[1]^(-1)
#
#def comp_syn_poly__(x, tildea:list, Ic:list, I:int, tildeg):
#  invs = [ lin_inv_mod( x-tildea[i], tildeg) for i in Ic if not i in I ]
#  if len(I)==1:
#    j = I[0]
#    prod_tmp = product(x-tildea[i] for i in Ic if i!=j).mod(tildeg)
#    return (x-tildea[j])*sum( prod_tmp*l_inv for l_inv in invs ) + prod_tmp
#  elif len(I)==2:
#    j1 = I[0]
#    j2 = I[1]
#    prod_tmp = product(x-tildea[i] for i in Ic if i!=j1 and i!=j2).mod(tildeg)
#    return (x-tildea[j1])*(x-tildea[j2]) * sum( (prod_tmp*l_inv).mod(tildeg) for l_inv in invs ) + (x-tildea[j1])*prod_tmp + (x-tildea[j2])*prod_tmp
#  else:
#    return comp_syn_poly(x, tildea, Ic)

"""
Computes the syndrome polynomial from partial support elements
Function needed for support completion if only a subset of support elements are known.
"""
#where index alpha[j] is expensive, i.e., only calculated later on!
def comp_syn_poly_(x, tildea:list, Ic:list, I:int):
  if len(I)==1:
    j = I[0]
    prod_tmp = product(x-tildea[i] for i in Ic if i!=j)
    return (x-tildea[j])*sum( prod_tmp.quo_rem(x-tildea[i])[0] for i in Ic if i!=j ) + prod_tmp
  elif len(I)==2:
    j1 = I[0]
    j2 = I[1]
    prod_tmp = product(x-tildea[i] for i in Ic if i!=j1 and i!=j2)
    return (x-tildea[j1])*(x-tildea[j2]) * sum( prod_tmp.quo_rem(x-tildea[i])[0] for i in Ic if i!=j1 and i!=j2 ) + (x-tildea[j1])*prod_tmp + (x-tildea[j2])*prod_tmp
  else:
    return comp_syn_poly(x, tildea, Ic)


def comp_ind_cw(G, rem_idxs:set, upd_idxs:set, guessed_idxs:set):
  return [ ([i], get_code_words(G, rem_idxs.union(upd_idxs) - set([i]) , set( [i]+list(guessed_idxs) ))[0]) for i in rem_idxs ]

def comp_goppa_cw(G, rem_idxs:set, upd_idxs:set, guessed_idxs:set):
  return get_code_words(G, rem_idxs.union(upd_idxs), guessed_idxs)[:5]

"""

Calculates the missing support elements using the computed Goppa polynomial
Completes a partial support tildea using tildeg.
This function is needed if the polynomial system contains only a subset of unknowns of the support set.
"""
def support_completion(tildea:list, tildeg, rem_idxs:set, upd_idxs:set, ind_cw):
  if len(rem_idxs)==0:
    return tildea

  print(f'start support completion')
  start_time = time.time()
  s_time = time.time()

  #construct rings
  F2ma = tildea[list(rem_idxs)[0]].parent()
  F2m = F2ma.base_ring()
  #F2ma = PolynomialRing(F2m, [tildea[i] for i in rem_idxs])
  F2max.<x> = PolynomialRing(F2ma, 'x')

  ##find ai's two-at-a-time -- NOT USED ANYMORE!
  #rem_idxs = list(rem_idxs)
  #compute polys for all c in cws and sums of consecutive pairs
  #cws = ind_cw + [ (ind_cw[i][0]+ind_cw[i-1][0], ind_cw[i][1]+ind_cw[i-1][1]) for i in range(1,len(rem_idxs)) ]
  cws = ind_cw
  cws = [ (I, list( i for i in range(len(c)) if c[i].is_one() ) ) for I,c in cws ]
  print(f'  assembled codewords (took {time.time()-start_time:.2f}s)')
  s_time = time.time()
  syn_polys = [ comp_syn_poly_(x, tildea, Ic, I) for I,Ic in cws ]
  #syn_polys == [ comp_syn_poly(x, tildea, Ic) for I,Ic in cws ]
  print(f'  computed syndrome polys (took {time.time()-start_time:.2f}s)')
  s_time = time.time()
  #extract eqs
  eqs = [ f for s_poly in syn_polys for f in s_poly.mod(tildeg).list() ]
  #+ [ f for s_poly in syn_polys for f in s_poly.mod(tildeg**2).list() ]
  I = F2ma.ideal( eqs )
  I = F2ma.ideal( I.interreduced_basis() )
  #assemble solution -- note I now must be generated by linear polys!
  if any(f.degree()>1 for f in I.gens()):
    raise ValueError('could not find unique extension to given tildea and with given code-words!')

  sol = dict()
  for f in I.gens():
    sol[f.lt()] = f.constant_coefficient()

  print(f'  computed solutions (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  #assert( len(sols)==1 )

  for i in rem_idxs:
    tildea[i] = tildea[i].subs(sol).constant_coefficient()

  for i in upd_idxs:
    tildea[i] = tildea[i].subs(sol).constant_coefficient()

  print(f'end support completion (took {time.time()-start_time:.2f}s)')
  return tildea


"""
Calls the C-program 'comparePK' to check if (tildea, tildeg) generate the correct key
Return: 1 on success, 0 on failure

"""
def verify_alt_gen(n,tildea, tildeg, fname_pk) -> bool:
  #write tildea
  fd, tmp_fname_a = tempfile.mkstemp()
  os.close(fd)
  write_alpha(tildea, tmp_fname_a)
  #write tildeg
  fd,tmp_fname_g = tempfile.mkstemp()
  os.close(fd)
  write_poly(tildeg, tmp_fname_g)
  #comparePK requires output_sage_a and output_sage_g to be written in its directory!
  sh_comparePK = './comparePK_'+str(n)+' ' + fname_pk + ' ' + tmp_fname_g + ' ' + tmp_fname_a + '> /dev/null 2>&1'
  output = os.system(sh_comparePK)
  return output==0
  #return output[-1] == 'SUCCESS: PUBLIC KEYS IDENTICAL '

"""
This function runs in parallel on different cores for all solutions found in the polynomial system of equations
It computes the Goppa polynomial, calculates (if needed) missing support elements of the support set
and checks if the calculated secret key pair is a valid alternative key
"""
@parallel(ncpus=NCORE)
def extend_support_candidate(n,t, tildea, rem_idxs,upd_idxs, ind_cw, goppa_cws):
  tildeg = goppaGCD(G, tildea, t, goppa_cws)
  if not tildeg.is_zero():
    try:
      tildeg = tildeg.nth_root(2)
      print('potential goppa-poly is a square, returning root of it!')
      #complete support
      tildea = support_completion(tildea, tildeg, rem_idxs, upd_idxs, ind_cw)
      s_time = time.time()
      if not verify_alt_gen(n,tildea, tildeg, fname_pk):
        print(f'alternative pair does NOT generate same code! -> continuing with next support candidate')
        return ([],tildeg.parent().zero())
      print(f'verified that alternative generators produce the same public key (took {time.time()-s_time:.2f}s)')
      return (tildea, tildeg)
    except ValueError:
      print('potential goppa-poly is NO square OR support completion failed! -> continuing with next support candidate')
      return ([],tildeg.parent().zero())

"""
This is the main function to compute an alternative secret key
Input: code length n, extension degree m, degree of Goppa polynomial t, generator matrix G, linear polynomials, quadratic polyinomials, filename of public key
Return: calculated support set tildea, calculated Goppa polynomial tildeg
"""
def compute_alternative_secret_key(n:int, m:int, t:int, G:list, lin_polys:list, quad_polys:list, fname_pk:str) -> tuple:
  #solve linear and quadratic polys
  lin_polys = solve_linear(n, lin_polys)
  S_L,rem_idxs,upd_idxs,guessed_idxs = solve_quad(n, m, lin_polys, quad_polys)

  s_time = time.time()
  #compute indicator-cw -- if rem_idxs!=[]
  ind_cw = comp_ind_cw(G, rem_idxs, upd_idxs, guessed_idxs)
  goppa_cw = comp_goppa_cw(G, rem_idxs, upd_idxs, guessed_idxs)
  print(f'computed \'indicator\' and \'goppaGCD\'-codewords (took {time.time()-s_time:.2f}s)')

  ###parallel processing of support candidate set
  total_processed = 0
  done = False
  while not done:
    #compute NCORE elemets from S_L
    tildeas = []
    #tmp = [next(S_L) for i in range(NCORE)]
    while len(tildeas)<NCORE and not done:
      try:
        tildeas.append( next(S_L) )
      except StopIteration:
        #if S_L is exhausted, stop while-loop after current iteration!
        done = True
    if NCORE==1:
      print(f'processing next support candidate, so far {total_processed} candidates have been checked.')
    else:
      print(f'processing next {len(tildeas)} support candidates, so far {total_processed} candidates have been checked.')
    #prepare inputs for parallel processing
    #outputs = list( extend_support_candidate( t,tildea,rem_idxs,upd_idxs,ind_cw ) for tildea in tildeas ) #serial implementation!
    inputs = [ (copy(n),copy(t),tildea,copy(rem_idxs),copy(upd_idxs),copy(ind_cw),copy(goppa_cw)) for tildea in tildeas ]
    outputs = list( extend_support_candidate( inputs ) )
    total_processed += len(tildeas)
    #check if one of outputs was successful!
    result = [o[1] for o in outputs if not o[1][1].is_zero()]
    if len(result)>0:
      return result[0]
  #for tildea in S_L:
  #  tildea, tildeg = extend_support_candidate(t,tildea,rem_idxs, ind_cw)
  #  if not tildeg.is_zero():
  #    return (tildea,tildeg)

  print('FAIL alternative goppa poly (with support) could not be found!')
  return ([],0)



if __name__=='__main__':
  #start timer
  s_time = time.time()
  start_time = time.time()
  
  #set up argparse
  parser = argparse.ArgumentParser(description='Solves fault equations and computes alternative secret key')
  parser.add_argument('-l', '--lin-polys', type=str, help='path to file containing linear polys; each line contains one linear polynomial represented by the corresponding indices, space-seperated; line "1 2 3 4" corresponds to polynomial "a[1]+a[2]+a[3]+a[4]".', required=True)
  parser.add_argument('-q', '--quad-polys', type=str, help='path to file containing quadratic polys; each line contains one quadratic polynomial represented by the corresponding indices, space-seperated; line "1 2 3 4" corresponds to polnyomial "a[1]*a[2] + a[3]*a[4]".', required=True)
  parser.add_argument('-pk', '--public-key', type=str, help='path to file containing public key (in hex-format), i.e., systematic part of PCM.', required=True)
  
  parser.add_argument('-a', '--alpha-out', type=str, help='path to file where found alternative support is written to; hex-representation, comma-separated.', default='output_sage_a.txt')
  parser.add_argument('-g', '--goppa-out', type=str, help='path to file where found alternative goppa-polynomial is written to; coeffs in hex-representation, comma-separated in increasing order (i.e. LC comes last)', default='output_sage_g.txt')
  
  parser.add_argument('-j', '--jobs', help='number of parallel threads for checking the support candidates', default=1, required=False, type=int)
  #parser.add_argument('-p', '--parallel', action='store_true', help='enable parallel checking of all support candidates')
  parser.add_argument('-n', help="code length", default=3488, required=False, type=int)
  parser.add_argument('-m', help="degree of field extension", default=12, required=False, type=int)
  parser.add_argument('-t', help="guaranteed error-correction capability", default=64, required=False, type=int)
  
  #parse args
  args = parser.parse_args()
  fname_l = args.lin_polys
  fname_q = args.quad_polys
  fname_pk = args.public_key

  fname_a = args.alpha_out
  fname_g = args.goppa_out

  #set actual number of cores -- should/must not exceed number of real processors!
  NCORE = min(NCORE, args.jobs)

  #check that comparePK binary is available and can be executed!
  if not os.access('./comparePK_'+str(args.n), os.X_OK):
    print('comparePK binary not found or not executable!')
    exit(1)
  
  #read linear and quadratic polys
  linear_polys = read_polys( fname_l )
  quadratic_polys = read_polys( fname_q )
  print(f'parsed fault eqs (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  G = code_gen_from_pk_hex( fname_pk, args.n, args.m, args.t )
  print(f'parsed public key and constructed generator matrix (took {time.time()-s_time:.2f}s)')
  s_time = time.time()
  
  ##find alternative secret key
  tildea,tildeg = compute_alternative_secret_key(args.n, args.m, args.t, G, linear_polys, quadratic_polys, fname_pk)
  if len(tildea)==0:
    print(f'alternative secret key could NOT be found (took {time.time()-s_time:.2f}s)')
    print(f'total time: {time.time()-start_time:.2f}s')
    exit(1)
  print(f'computed alternative secret key (took {time.time()-s_time:.2f}s)')
  print(f'writing alternative generators to \'{fname_a}\' and \'{fname_g}\'.')
  #write polys to file
  write_alpha(tildea, fname_a)
  write_poly(tildeg, fname_g)
  print(f'total time: {time.time()-start_time:.2f}s')



