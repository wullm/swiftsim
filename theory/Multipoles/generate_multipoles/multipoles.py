import numpy as np
import sys
import scipy.special

def factorial(x):
    if x == 0:
        return 1
    else:
        return x * factorial(x-1)

SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal(num):
    suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

# Get the order
order = int(sys.argv[1])

print "-------------------------------------------------"
print "Generating code for multipoles of order", order, "(only)."
print "-------------------------------------------------\n"

print "-------------------------------------------------"

print "Field tensor structure:"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
count = 0
tmp_str = '  float'
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                tmp_str = tmp_str + " F_%d%d%d,"%(i,j,k)
                count += 1
                if count == 3:
                    tmp_str = tmp_str[:-1]
                    print tmp_str + ';'
                    count = 0
                    tmp_str = '  float'
if order > 0:
    print "#endif"

print ""

print "-------------------------------------------------"
print "Multipole structure:"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
count = 0
tmp_str = '  float'
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                tmp_str = tmp_str + " M_%d%d%d,"%(i,j,k)
                count += 1
                if count == 3:
                    tmp_str = tmp_str[:-1]
                    print tmp_str + ';'
                    count = 0
                    tmp_str = '  float'
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_field_tensors_add():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  la->F_%d%d%d += lb->F_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_multipole_add():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  ma->M_%d%d%d += mb->M_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_multipole_equal():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

# Create all the terms relevent for this order
print "  /* Manhattan Norm of %s order terms */"%ordinal(order)
print "  const float order%d_norm =\n"%order,
first = True
count = 0
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                count += 1
                if first:
                    first = False
                    print "     ",
                else:
                    print "    +",
                print "fabsf(ma->M_%d%d%d)"%(i,j,k),
                if count == scipy.special.binom(order+2, 2):
                    print "+ fabsf(mb->M_%d%d%d);"%(i,j,k)
                else:
                    print "+ fabsf(mb->M_%d%d%d)"%(i,j,k)
print '\n'
print "  /* Compare %s order terms above 1%% of norm */"%ordinal(order)
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  if (fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > 0.01f * order%d_norm &&"%(i,j,k,i,j,k,order)
                print "      fabsf(ma->M_%d%d%d - mb->M_%d%d%d) / fabsf(ma->M_%d%d%d + mb->M_%d%d%d) > tolerance) {"%(i,j,k,i,j,k,i,j,k,i,j,k)
                print "    message(\"M_%d%d%d term different\");"%(i,j,k)
                print "    return 0;"
                print "  }"

if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"
print "gravity_P2M(): (init)"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)
# Create all the terms relevent for this order
count = 0
tmp_str = '  double'
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                tmp_str = tmp_str + " M_%d%d%d = 0.,"%(i,j,k)
                count += 1
                if count == 3:
                    tmp_str = tmp_str[:-1]
                    print tmp_str + ';'
                    count = 0
                    tmp_str = '  double'
if order > 0:
    print "#endif"
print ""
print "-------------------------------------------------"
print "gravity_P2M(): (loop)"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                if order % 2 == 0:
                    print "    M_%d%d%d += m * X_%d%d%d(dx);"%(i,j,k,i,j,k)
                else:
                    print "    M_%d%d%d += -m * X_%d%d%d(dx);"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"
print ""
print "-------------------------------------------------"
    
print "gravity_P2M(): (storing)"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "  /* %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  multi->m_pole.M_%d%d%d = M_%d%d%d;"%(i,j,k,i,j,k)
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_M2M():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d"%(order-1)

print "  /* Shift %s order terms */"%ordinal(order)
    
# Create all the terms relevent for this order
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                final_count = 0
                for iiii in range(2):
                    if iiii == 1: print "  m_a->M_%d%d%d =\n    m_b->M_%d%d%d"%(i,j,k,i,j,k),
                    count = 0
                    first = True
                    for ii in range(order+1):
                        for jj in range(order+1):
                            for kk in range(order+1):

                                if not(ii == 0 and jj == 0 and kk == 0):
                                    for iii in range(order+1):
                                        for jjj in range(order+1):
                                            for kkk in range(order+1):
                                                if ii+iii == i and jj+jjj == j and kk+kkk == k:
                                                    count += 1
                                                    if iiii == 0: final_count += 1
                                                    if count % 2 == 0:
                                                        if iiii == 1:
                                                            if count == final_count:
                                                                print "+ X_%d%d%d(dx) * m_b->M_%d%d%d;"%(ii, jj, kk, iii, jjj, kkk)
                                                            else:
                                                                print "+ X_%d%d%d(dx) * m_b->M_%d%d%d"%(ii, jj, kk, iii, jjj, kkk)
                                                    else:
                                                        if first:
                                                            if iiii == 1:
                                                                if count == final_count:
                                                                    print "+ X_%d%d%d(dx) * m_b->M_%d%d%d;"%(ii, jj, kk, iii, jjj, kkk)
                                                                else:
                                                                    print "+ X_%d%d%d(dx) * m_b->M_%d%d%d"%(ii, jj, kk, iii, jjj, kkk),
                                                            first = False
                                                        else:
                                                            if iiii == 1:
                                                                if count == final_count:
                                                                    print "    + X_%d%d%d(dx) * m_b->M_%d%d%d;"%(ii, jj, kk, iii, jjj, kkk)
                                                                else:
                                                                    print "    + X_%d%d%d(dx) * m_b->M_%d%d%d"%(ii, jj, kk, iii, jjj, kkk),

if order > 0:
    print "#endif"
    
print ""
print "-------------------------------------------------"

print "gravity_M2L():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

assert order > 0, 'M2L not right for order = 0'
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  const float M_%d%d%d = m_a->M_%d%d%d;"%(i,j,k,i,j,k)
print ""

for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                print "  const float D_%d%d%d = pot->D_%d%d%d;"%(i,j,k,i,j,k)
print ""

# Loop over LHS order
for l in range(order + 1):
    print "  /* Compute %s order field tensor terms (addition to rank %d) */"%(ordinal(order), l)

    for i in range(l+1):
        for j in range(l+1):
            for k in range(l+1):
                if i + j + k == l:
                    print "  l_b->F_%d%d%d +="%(i,j,k),
                    final_count = 0
                    for iiii in range(2):
                        count = 0
                        first = True
                        for ii in range(order+1):
                            for jj in range(order+1):
                                for kk in range(order+1):
                                    if ii + jj + kk  == order - l:
                                        count += 1
                                        if iiii == 0: final_count += 1
                                        if first:
                                            first = False
                                        else:
                                            if iiii == 1:
                                                if (count + 2)% 3 == 0:
                                                    print "    +",
                                                else:
                                                    print "+",
                                        if iiii == 1:
                                            if count == final_count:
                                                print "M_%d%d%d * D_%d%d%d;"%(ii,jj,kk,i+ii,j+jj,k+kk)
                                            else:
                                                if count % 3 == 0:
                                                    print "M_%d%d%d * D_%d%d%d"%(ii,jj,kk,i+ii,j+jj,k+kk)
                                                else:
                                                    print "M_%d%d%d * D_%d%d%d"%(ii,jj,kk,i+ii,j+jj,k+kk),
    print ""
    
if order > 0:
    print "#endif"
print ""
print "-------------------------------------------------"

print "gravity_L2L():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

# Loop over LHS order
for l in range(order + 1):
    print "  /* Shift %s order field tensor terms (addition to rank %d) */"%(ordinal(order), l)

    for i in range(l+1):
        for j in range(l+1):
            for k in range(l+1):
                if i + j + k == l:
                    print "  la->F_%d%d%d += "%(i,j,k),
                    final_count = 0
                    for iiii in range(2):
                        count = 0
                        first = True
                        for ii in range(order+1):
                            for jj in range(order+1):
                                for kk in range(order+1):
                                    if ii + jj + kk  == order - l:
                                        count += 1
                                        if iiii == 0: final_count += 1
                                        if first:
                                            first = False
                                        else:
                                            if iiii == 1 and count % 2 != 0: print "    +",
                                        if iiii == 1:
                                            if count == final_count:
                                                print "X_%d%d%d(dx) * lb->F_%d%d%d;"%(ii,jj,kk,i+ii,j+jj,k+kk)
                                            else:
                                                if count % 2 == 0:
                                                    print "X_%d%d%d(dx) * lb->F_%d%d%d"%(ii,jj,kk,i+ii,j+jj,k+kk)
                                                else:
                                                    print "X_%d%d%d(dx) * lb->F_%d%d%d"%(ii,jj,kk,i+ii,j+jj,k+kk),
    print ""
if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

print "gravity_L2P():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)

    print "  /* %s order contributions */"%(ordinal(order))

    for r in range(3):
        print "  a_grav[%d] +="%(r),

        supercount = 0
        for iiii in range(2):
            first = True
            count = 0
            for i in range(order + 1):
                for j in range(order + 1):
                    for k in range(order + 1):
                        if i + j + k == order-1:
                            count += 1
                            if iiii == 0: supercount += 1
                            if first:
                                first = False
                            else:
                                if count %2 == 0:
                                    if iiii == 1: print "+",
                                else:
                                    if iiii == 1: print "    +",
                            if r == 0:
                                ii = i+1
                                jj = j
                                kk = k
                            if r == 1:
                                ii = i
                                jj = j+1
                                kk = k
                            if r == 2:
                                ii = i
                                jj = j
                                kk = k+1
                            if iiii == 1:
                                if count % 2 == 0:
                                    if count == supercount:
                                        print "X_%d%d%d(dx) * lb->F_%d%d%d;"%(i,j,k,ii,jj,kk)
                                    else:
                                        print "X_%d%d%d(dx) * lb->F_%d%d%d"%(i,j,k,ii,jj,kk)
                                else:
                                    if count == supercount:
                                        print "X_%d%d%d(dx) * lb->F_%d%d%d;"%(i,j,k,ii,jj,kk)
                                    else:
                                        print "X_%d%d%d(dx) * lb->F_%d%d%d"%(i,j,k,ii,jj,kk),

    print ""

# Contribution to particle potential.
print "  pot -=",
count = 0
first = True
for i in range(order+1):
    for j in range(order+1):
        for k in range(order+1):
            if i + j + k == order:
                count += 1
                if first:
                    first = False
                else:
                    if count % 2 != 0:
                        print "    +",
                    else:
                        print "+",
                if count % 2 == 0:
                    if count == scipy.special.binom(order+2, 2):
                        print "X_%d%d%d(dx) * lb->F_%d%d%d;"%(i,j,k,i,j,k)
                    else:
                        print "X_%d%d%d(dx) * lb->F_%d%d%d"%(i,j,k,i,j,k)
                else:
                    if count == scipy.special.binom(order+2, 2):
                        print "X_%d%d%d(dx) * lb->F_%d%d%d;"%(i,j,k,i,j,k)
                    else:
                        print "X_%d%d%d(dx) * lb->F_%d%d%d"%(i,j,k,i,j,k),
print ""
if order > 0:
    print "#endif"
print ""
print "-------------------------------------------------"

print "gravity_M2P():"
print "-------------------------------------------------\n"

if order > 0:
    print "#if SELF_GRAVITY_MULTIPOLE_ORDER > %d\n"%(order-1)
    
print "  /* %s order contributions */"%(ordinal(order))

    
for r in range(4):
    if r == 0:
        print "  *f_x =",
    if r == 1:
        print "  *f_y =",
    if r == 2:
        print "  *f_z =",
    if r == 3:
        print "  *pot =",

    supercount = 0
    for iiii in range(2):     
        first = True
        count = 0
        for i in range(order+1):
            for j in range(order+1):
                for k in range(order+1):
                    if i + j + k == order:
                        count += 1
                        if iiii == 0: supercount += 1
                        if first:
                            first = False
                        else:
                            if count %2 == 0:
                                if iiii == 1: print "+",
                            else:
                                if iiii == 1: print "    +",
                        if r == 0:
                            ii = i+1
                            jj = j
                            kk = k
                        if r == 1:
                            ii = i
                            jj = j+1
                            kk = k
                        if r == 2:
                            ii = i
                            jj = j
                            kk = k+1
                        if r == 3:
                            ii = i
                            jj = j
                            kk = k
                        if iiii == 1:
                            if count % 2 == 0:
                                if count == supercount:
                                    print "m->M_%d%d%d * d.D_%d%d%d;"%(i,j,k,ii,jj,kk)
                                else:
                                    print "m->M_%d%d%d * d.D_%d%d%d"%(i,j,k,ii,jj,kk) 
                            else:
                                if count == supercount:
                                    print "m->M_%d%d%d * d.D_%d%d%d;"%(i,j,k,ii,jj,kk)
                                else:
                                    print "m->M_%d%d%d * d.D_%d%d%d"%(i,j,k,ii,jj,kk),
                    
print ""

if order > 0:
    print "#endif"

print ""
print "-------------------------------------------------"

