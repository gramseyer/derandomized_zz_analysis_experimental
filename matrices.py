from sage.all import *

import numpy
from numpy import linalg


def add_lambda(x, y, generator):
   if (x ^ generator == y):
      return 1 
   else:
      return 0

def small_matrix(d1, generator):
    return numpy.array(matrix(RR, d1, d1, lambda i, j: add_lambda(i,j,generator)))

def big_matrix(n, d1, generator):
  return numpy.kron(numpy.eye(n), small_matrix(d1, generator))

def parity(z):
  counter = 0
  for q in z.binary():
    if (q == '1'):
      counter = counter + 1
  return counter % 2

def bin_inner_prod(a, b):
  return a & b 

def bin_string(x, y,  n):
   string = ''
   for i in range(0, n):
       xi = x**i
       coeffs = xi.polynomial().coefficients(sparse=False)
       coeff_str = '+0b' + ''.join(map(str, coeffs[::-1])).zfill(n)
       result = Integer(coeff_str)
       bitstr = bin_inner_prod(y, result)
       bit = parity(bitstr)
       string = string + str(bit)
   #print string
   return string

def binary_map(y):
  return {
    '00' : '00',
    '01' : '10',
    '10' : '11',
    '11' : '01',
  }[y]

def ternary_map(y):
  return {
    '000' : '000',
    '001' : '010',
    '010' : '100',
    '011' : '110',
    '100' : '011',
    '101' : '001',
    '110' : '111',
    '111' : '101',
  }[y]

def gen_map(y, m):
  if m == 0:
    return '' 
  if m % 2 == 1:
    return gen_map(y // 8, m-3) + ternary_map ("{0:b}".format(y % 8).zfill(3))
  return gen_map(y//4, m-2) + binary_map("{0:b}".format(y % 4).zfill(2))

def generator_pair(x, y, m, n):
  str1 = bin_string(x, y, n)
  str2 = bin_string(x, Integer('+0b' + gen_map(y, m)), n)
  return (str1, str2)

def gen_graph(n, d1):
  g = graphs.RandomRegular(d1, n)
  edge_colorings = graph_coloring.edge_coloring(g)
  if not (len(edge_colorings) == d1):
    return None
  return gen_perm_mat(n, d1, g, edge_colorings)

def gen_perm_mat(n, d1, g, edge_colorings):
  perm_matrix = matrix(RR, n*d1)
  for i in range(0, len(edge_colorings)):
    edge_coloring = edge_colorings[i]
    for (a,b) in edge_coloring:
      perm_matrix[a*d1 + i, b*d1 + i] = 1
      perm_matrix[b*d1 + i, a*d1 + i] = 1
  return perm_matrix

def get_single_lin_operator(gen1, gen2, num_vertex_clouds, d1, perm_mat):
  return numpy.matmul(numpy.matmul(big_matrix(num_vertex_clouds, d1, gen2), perm_mat), big_matrix(num_vertex_clouds, d1, gen1))

def get_entire_lin_operator(generators, num_vertex_clouds, d1, perm_mat):
  mat = numpy.zeros((num_vertex_clouds * d1, num_vertex_clouds * d1)) #matrix(RR, num_vertex_clouds * d1)
  count = 0
  for (gen1, gen2) in generators:
    mat = mat + get_single_lin_operator(gen1, gen2, num_vertex_clouds, d1, perm_mat)
    count = count + 1
    if (count % (len(generators)/10) == 0):
      print count
  return mat / len(generators)

def gen_all_generators(m, n_gen):
  gen_pairs = []
  for x in GF(2**m):
    for y in range(0, 2**m):
      (g1, g2) = generator_pair(x, y, m, n_gen)
      #print (g1, g2)
      gen_pairs = gen_pairs + [(Integer('+0b' + g1), Integer('+0b' + g2))]
  return gen_pairs

def get_individual_matrix(d1, generators):
  mat = matrix(QQ, d1)
  for gen in generators:
    mat = mat + small_matrix(d1, gen)
  mat = mat / len(generators)
  return mat

def get_correlated_matrix(d1, generators):
  mat = matrix(QQ, d1)
  for (gen, gen2) in generators:
    mat = mat + (small_matrix(d1, gen2) * small_matrix(d1, gen))
  mat = mat / len(generators)
  return mat

def test_bias(test, generators):
  zero = 0.0
  one = 0.0
  for g1 in generators:
     if (parity(bin_inner_prod(g1, test)) == 1):
        one = one + 1
     else:
        zero = zero + 1
  return (one - zero)/ (zero + one)

def fst(x):
  return x[0]

def snd(x):
  print "foo"
  return x[1]

def fstsnd(x):
  return (x[0]).__xor__(x[1])

def test_all_bias(generators, max_test, which=fst):
  max_bias = 0
  for i in range(1, max_test):
    test_result = abs(test_bias(i, (map(which, generators))))
    if (test_result > max_bias):
      max_bias = test_result
  print ("resulting bias is " + str(max_bias))
  return max_bias

def get_random_block_perp(d1):
  acc = Rational(0)
  vec = []
  for i in range(0, d1-1):
    rand = RR.random_element(-10, 10)
    acc = acc + rand
    vec = vec + [rand]
  vec = vec + [-acc]
  return vec


def get_random_alphaperp_vector(num_vertex_clouds, d1):
  vec = []
  for i in range(0, num_vertex_clouds):
    vec = vec + get_random_block_perp(d1)
  return vector(vec)

def run_random_ap_test(matrix, num_vertex_clouds, d1, lambda2):
  ap = get_random_alphaperp_vector(num_vertex_clouds, d1)
  result = matrix * ap
  ratio = (result.norm())/(ap.norm())
  #ratio = (abs(result.dot_product(ap))) / (ap.dot_product(ap))
  if (ratio > lambda2):
    print (ratio.n())
  return (ratio, ap)

def run_endless_tests(matrix, num_vertex_clouds, d1, lambda2, test_max=10000):
  current_ratio = lambda2
  test_count = 0
  bad_vecs = []
  print ("performing " + str(test_max) + " tests")
  while(True):
    test_count = test_count + 1
    (ratio, vec) = run_random_ap_test(matrix, num_vertex_clouds, d1, current_ratio)
    if ratio > current_ratio:
      current_ratio = ratio
      bad_vecs = bad_vecs + [vec]
      print (ratio.n(), vec)
    if test_count % (test_max/5) == 0:
      print str ((100*test_count/test_max)) + "%"
    if (test_count > test_max):
      return bad_vecs

def gen_good_matrix(num_clouds, d1, epsilon):
  perm_mat, g, g_eigen = gen_uncolorable_graph(num_clouds, d1)
  return test_good_matrix(num_clouds, d1, epsilon, perm_mat, g)  

def test_good_matrix(num_clouds, d1, epsilon, perm_mat, g):
  gen_n = log(d1, 2)
  if (gen_n - ceil(gen_n) != 0):
    print "d1 must be power of 2"
    return None
  gen_m = ceil (log(gen_n/epsilon, 2))
  print(gen_m, gen_n)
  generators = gen_all_generators(gen_m, gen_n)
  d2 = len(generators)
  print ("d2 = " + str(d2))
  bias = test_all_bias(generators, d1)
  #(perm_mat, g, g_eigen) = gen_uncolorable_graph(num_clouds, d1)
  #print ("made perm mat")
  total_matrix = get_entire_lin_operator(generators, num_clouds, d1, perm_mat)
  print("Generated matrix with lambda2 = " + str(bias))
  return (total_matrix, bias, g)

def run_random_test(num_clouds, d1, epsilon):
  (total_matrix, bias, g) = gen_good_matrix(num_clouds, d1, epsilon)
  bad_vecs = run_endless_tests(total_matrix, num_clouds, d1, bias)
  return (total_matrix, bad_vecs, g)

def run_until_failure(num_clouds, d1, epsilon):
  while True:
    results = run_random_test(num_clouds, d1, epsilon)
    if not (len(results[1]) == 0):
      return results

def second_eigenval(total_matrix):
  eigenvalues = total_matrix.eigenvalues()
  tosort = map(abs, eigenvalues)
  list.sort(tosort, reverse = True)
  return tosort[1]

def normed_eigenvals(total_matrix):
  eigenvalues = total_matrix.eigenvalues()
  tosort = map (abs, eigenvalues)
  list.sort(tosort, reverse=True)
  return tosort

def gen_colorable_graph(num_clouds, d1):
  g = None
  ec = None
  while True:
    g = graphs.RandomRegular(d1, num_clouds)
    ec = graph_coloring.edge_coloring(g)
    if len(ec) == d1:
      break
    print "picking failed, trying again"
  p = gen_perm_mat(num_clouds, d1, g, ec)
  return (p, g, 1) #1 is placeholder

def gen_uncolorable_graph(num_clouds, d1):
  g = graphs.RandomRegular(d1, num_clouds)
  print "picked a random regular graph"
  return gen_uncolorable_graph_perm_mat(num_clouds, d1, g)

def gen_uncolorable_graph_perm_mat(num_clouds, d1, g):
  edge_list = g.to_dictionary()
  mat = numpy.zeros((num_clouds * d1, num_clouds * d1)) #matrix(RR, num_clouds * d1)
  for (cloud, edges) in edge_list.iteritems():
    for i in range (0, d1):
      list.sort(edges)
      target = edges[i]
      incoming_edge = edge_list[target].index(cloud)
      mat[cloud*d1 + i, target*d1 + incoming_edge] = 1
  eigen = 0 #second_eigenval(g.am()/d1)
  return (mat, g, eigen)


#index is 1 through length-1 inclusive
def gen_ortho_ap_basis_vec(length, index):
  return ([1/(sqrt(index + index*index))] * index + [-index/sqrt(index + index * index)]) + ([0] * (length - index - 1))

def gen_ortho_ap_basis_submatrix(length):
  mat = numpy.zeros((length, length-1))
  for j in range(0, length-1):
    vec = gen_ortho_ap_basis_vec(length, j+1)
    for i in range(0, length):
      mat[(i, j)] = vec[i]
  return mat

def gen_ortho_ap_basis_matrix(num_clouds, d_1):
  return numpy.kron(numpy.eye(num_clouds), gen_ortho_ap_basis_submatrix(d_1)) 

def get_random_numpy_array(num_clouds, d1, bias):
  (matrix, bias, g) = gen_good_matrix(num_clouds, d1, bias)
  proj = gen_ortho_ap_basis_matrix(num_clouds, d1)
  return (numpy.matmul(numpy.matmul(numpy.transpose(proj), matrix), proj), bias, g)

def test_numpy_array(num_clouds, d1, matrix):
  proj = gen_ortho_ap_basis_matrix(num_clouds, d1)
  return numpy.matmul(numpy.matmul(numpy.transpose(proj), matrix), proj)

def get_eigenvals_vecs(mat):
  e_vals, e_vecs = linalg.eig(mat)
  abs_vec = numpy.vectorize(abs)
  e_vals_abs = abs_vec(e_vals)
  idx = e_vals_abs.argsort()[::-1]
  return (e_vals_abs[idx], e_vecs[:,idx])

def get_largest_eig(mat):
  vals, vecs = get_eigenvals_vecs(mat)
  return (vals[0], vecs[0])

def run_many_eig_tests(num_clouds, d1, epsilon):
  log = []
  graphs = []
  goal_bias = 0
  while True:
  #for i in range (0, 100): 
    (mat, bias, g) = get_random_numpy_array(num_clouds, d1, epsilon)
    actual_bias = get_largest_eig(mat)[0]
    log = [actual_bias]
    graphs = [g]
    goal_bias = bias
    print actual_bias
    if actual_bias < bias + 0.00001:
      break
  return (log, graphs, goal_bias)

def filter_good(mats, graphs):
  output = []
  for i in range (0, len(mats)):
    if mats[i] < 0.126:
      output = output + [(mats[i], graphs[i])]
  return output

def filter_bad(mats, graphs):
  output = []
  for i in range(0, len(mats)):
    if mats[i] > 0.126:
      output = output + [graphs[i]]
  return output

def gen_super_graph(num_clouds, d1):
  if not (log(d1, 2) == ceil (log (d1, 2))):
    print "d1 must be power of two"
    return None
  half = d1/2
  edges = {}
  for i in range (0, num_clouds/half):
    for j in range(0, half):
      edges[i*half + j] = map(lambda x : mod(x, num_clouds), list(xrange((i+1)*half, (i+2)*half)))
  return Graph(edges)
  
def super_test(num_clouds, d1, bias):
  g = gen_super_graph(num_clouds, d1)
  (p_mat, _, _) = gen_uncolorable_graph_perm_mat(num_clouds, d1, g)
  (total_mat, goal_bias, _) = test_good_matrix(num_clouds, d1, bias, p_mat, g)
  return get_largest_eig(test_numpy_array(num_clouds, d1, total_mat))[0]

def general_test(num_clouds, d1, bias, g):
  p_mat, _, _ = gen_uncolorable_graph_perm_mat(num_clouds, d1, g)
  (total_mat, goal_bias, _) = test_good_matrix(num_clouds, d1, bias, p_mat, g)
  return get_largest_eig(test_numpy_array(num_clouds, d1, total_mat))

def gen_div_list(factor, list_items):
  output = []
  list.sort(list_items)
  for i in range(0, len(list_items)):
    output = output + list(xrange(factor * list_items[i], factor * list_items[i] + factor))
  return output

def gen_div_graph(num_clouds, d1, factor):
  g = graphs.RandomRegular(Integer(d1/factor), num_clouds/factor)
  return dup_graph(factor, g)

def dup_graph(factor, g):
  d = {k : gen_div_list(factor, sorted(v)) for k, v in g.to_dictionary().items()}
  results = []
  for (k, v) in d.iteritems():
    for i in range(0, factor):
      results = results + [(factor *k + i, v)]
  return Graph((dict(results)))


