#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time

import numpy as np
import pyopencl as cl
from scipy.misc import imread
from matplotlib import pyplot as plt
from envmap.environmentmap import EnvironmentMap


def associatedLegendrePolynomialValues(mm, n, m, x):
    """
    mm: number of points
    n: order, 0 <= n
    m: mode (left-right), 0 <= m <= n
    x: function (of size mm, minimum)
    References:
    https://people.sc.fsu.edu/~jburkardt/m_src/legendre_polynomial/pm_polynomial_value.m
    https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    """
    cx = np.zeros((mm, n+1))

    if m <= n:
        cx[:mm, m] = 1.

        fact = 1.
        for _ in range(m):
            cx[:mm, m] = - cx[:mm, m] * fact * np.sqrt( 1. - x[:mm]**2 )
            fact = fact + 2.

    if m + 1 <= n:
        cx[:mm, m+1] = ( 2*m + 1 ) * x[:mm] * cx[:mm, m]

    for j in range(m + 2, n + 1):
        cx[:mm, j] = ( ( ( 2 * j     - 1 ) * x[:mm] * cx[:mm, j - 1]
                       + (   - j - m + 1 ) * cx[:mm, j - 2] )
                       / (     j - m     ) )

    return cx

def spharm(n, m, u, v):
    """
    https://imdoingitwrong.wordpress.com/2011/04/14/spherical-harmonics-wtf/
    """
    N = np.sqrt( (2*n + 1) / (4 * np.pi) * (
        np.math.factorial(n - np.abs(m)) / np.math.factorial(n + np.abs(m))))
    if m != 0:
        N *= np.sqrt(2)
    Phi = np.exp(1j*np.abs(m)*v) * (-1)**((m<0)*m)
    P = associatedLegendrePolynomialValues(u.shape[0], n, np.abs(m), u)[:,-1]
    return N * np.outer(P, Phi)


def sphericalHarmonics(order):
    e = EnvironmentMap(50, 'LatLong')
    sa = e.solidAngles()

    cols = ['Mode = {}'.format(col) for col in range(-order, order+1)]
    rows = ['Order = {}'.format(row) for row in range(order+1)]

    u = np.cos(np.linspace(0, np.pi, 50))
    v = np.linspace(0, 2*np.pi, 100)
    ret = []
    for n in range(order):
        for m in range(-n, n + 1):
            s = spharm(n, m, u, v)
            ret.append(s * sa)
            plt.figure(0)
            plt.subplot(order, order*2-1, n*(order*2-1)+m+order)
            plt.imshow(ret[-1].real)
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(False)
            #plt.title('n: {}, m: {}'.format(n, m))
            plt.figure(1)
            plt.subplot(order, order*2-1, n*(order*2-1)+m+order)
            plt.imshow(ret[-1].imag)
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(False)
            #plt.title('n: {}, m: {}'.format(n, m))
    plt.figure(0)
    plt.tight_layout()
    #plt.gcf().suptitle("order (m)", fontsize=14)
    #plt.gcf().text(0.02, 0.5, 'degree (n)',
    #               horizontalalignment='center',
    #               verticalalignment='top', fontsize=14,
    #               rotation='vertical')
    fig = plt.gcf()
    fig.set_size_inches(18.5, 5.2)
    fig.savefig('real.png', dpi=100, bbox_inches='tight')

    plt.figure(1)
    plt.tight_layout()
    #plt.gcf().suptitle("order (m)", fontsize=14)
    #plt.gcf().text(0.02, 0.5, 'degree (n)',
    #               horizontalalignment='center',
    #               verticalalignment='top', fontsize=14,
    #               rotation='vertical')
    fig = plt.gcf()
    fig.set_size_inches(18.5, 5.2)
    fig.savefig('imag.png', dpi=100, bbox_inches='tight')
    
    #plt.show()
    return ret


def spharmOCL(im):
    assert im.shape[0]*2 == im.shape[1]

    e = EnvironmentMap(im.shape[0], 'LatLong')
    sa = e.solidAngles().astype('float32')
    #plt.imshow(sa)
    #plt.show()
    #return

    depth = max(im.shape)
    res_r = np.empty(depth**2, dtype='float32')
    res_i = np.empty(depth**2, dtype='float32')

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    mf = cl.mem_flags
    sa_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=sa)
    im_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=im)

    res_r_g = cl.Buffer(ctx, mf.READ_WRITE, res_r.nbytes)
    res_i_g = cl.Buffer(ctx, mf.READ_WRITE, res_i.nbytes)

    prg = cl.Program(ctx, """
void atomic_add_global(volatile global float *source, const float operand);

inline float divfact(uint a, uint b) {
    float r = 1.0f;
    for (uint i = a + 1; i <= b; ++i) {
        r *= i;
    }
    return 1.0f/r;
}

void atomic_add_global(volatile global float *source, const float operand) {
    union {
        unsigned int intVal;
        float floatVal;
    } newVal;
    union {
        unsigned int intVal;
        float floatVal;
    } prevVal;
 
    do {
        prevVal.floatVal = *source;
        newVal.floatVal = prevVal.floatVal + operand;
    } while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

__kernel void spharm(__global const uchar *im,
                  const uint height,
                  const uint width,
                  const uint depth,
                  __global const float *sa,
                  __global float *res_r,
                  __global float *res_i) {
    const uint gidy = get_global_id(0);
    const uint gidx = get_global_id(1);
    const float this_sa = sa[uint(gidy*width)+gidx];
    const float px = convert_float_rte(im[uint(gidy*width)+gidx]);

    if (gidx >= width || gidy >= height) {
        return;
    }
    if (gidx != 1 || gidy != 1) return;

    float v = (float)gidx / ((float)width - 1.0f) * 2.0f * M_PI;
    float u = cos((float)gidy / ((float)height - 1.0f) * M_PI);

    //uint nmax = width/2;
    uint nmax = depth;
    float N, Phi_r_p, Phi_i_p, Phi_r_n, Phi_i_n, P, P_1, P_2;
    float fact;
    int s_index = 0;

    for (uint m = 0; m < nmax; ++m) {
        Phi_r_p = cos((float)m*v);
        Phi_i_p = sin((float)m*v);
        Phi_r_n = cos(-(float)m*v);
        Phi_i_n = sin(-(float)m*v);

        P_2 = 1.0f;
        fact = 1.0f;
        for (uint i = 0; i < m; ++i) {
            P_2 = -P_2 * fact * sqrt( 1.0f - (u * u) );
            fact = fact + 2.0f;
        }
            if (isnan(P_2) || isinf(P_2)) {
                printf("Here 1! nmax: %u\\n", nmax);
                printf("P-2: %f, P-1: %f, P: %f\\n", P_2, P_1, P);
                return;
            }

        P_1 = ( 2*m + 1 ) * u * P_2;

            if (isnan(P_1) || isinf(P_1)) {
                printf("Here 2! nmax: %u\\n", nmax);
                printf("P-2: %f, P-1: %f, P: %f\\n", P_2, P_1, P);
                return;
            }

        for (uint n = m; n < nmax; ++n) {
            // sqrt(2) ? Is it a normalization term?
            N = sqrt( (float)(2 >> (m == 0)) * (2*n + 1) / (4.0f * M_PI) * (divfact(n - m, n + m)));

            P = ( ( ( 2.0f * ((float)n + 2.0f)            - 1.0f ) * u * P_1
                  + (      - ((float)n + 2.0f) - (float)m + 1.0f ) * P_2 )
                  / (        ((float)n + 2.0f) - (float)m        ) );

            //printf("N: %f, Phi_r_p: %f, Phi_i_p: %f, P_2: %f, sa: %f, px: %f, divfact: %f \\n", N, Phi_r_p, Phi_i_p, P_2, this_sa, px, divfact(n-m, n+m));
            //printf("%u: %f, ", s_index, res_r[s_index]);
            //printf("P-2: %f, P-1: %f, P: %f\\n", P_2, P_1, P);
            //printf("-> %f\\n", (2.0f*((float)n+2.0f) - 1.0f)*u);
            //printf("u: %f, n: %f, part1: %f, part2: %f, num: %f, den: %f\\n", u, (float)n, (2.0f*((float)n+2.0f)-1.0f)*u*P_1, (-(n+2.0f) -m + 1)*P_2, (2.0f*((float)n+2.0f)-1.0f)*u*P_1 + (-((float)n+2.0f) -(float)m + 1.0f)*P_2, ((float)n+2.0f)-(float)m);
            if (isnan(P) || isinf(P) || isnan(N*Phi_r_n*P_2*this_sa*px)) {
                printf("nmax: %u\\n", nmax);
                printf("P-2: %f, P-1: %f, P: %f\\n", P_2, P_1, P);
                printf("u: %f, n: %f, m: %f, part1: %f, part2: %f, num: %f, den: %f\\n", u, (float)n, (float)m, (2.0f*((float)n+2.0f)-1.0f)*u*P_1, (-(n+2.0f) -m + 1)*P_2, (2.0f*((float)n+2.0f)-1.0f)*u*P_1 + (-((float)n+2.0f) -(float)m + 1.0f)*P_2, ((float)n+2.0f)-(float)m);
                return;
            }
            
            atomic_add_global(res_r+s_index, N*Phi_r_n*P_2*this_sa*px);
            //printf("%f\\n", res_r[s_index]);
            atomic_add_global(res_i+s_index, N*Phi_i_n*P_2*this_sa*px);
            if (m != 0) {
                atomic_add_global(res_r+s_index+1, N*Phi_r_p*P_2*this_sa*px);
                atomic_add_global(res_i+s_index+1, N*Phi_i_p*P_2*this_sa*px);
            }
            s_index += (2 >> (m == 0)); // +1 if m == 0, +2 otherwise

            P_2 = P_1;
            P_1 = P;
        }
    }
    //printf("Pixel finished!\\n");
}
    """).build()

    global_shape = tuple(map(int, np.ceil(np.array(im.shape)/64.0)*64))
    prg.spharm(queue, global_shape, (1, 1),
               im_g,
               np.uint32(im.shape[0]),
               np.uint32(im.shape[1]),
               np.uint32(10),
               sa_g,
               res_r_g,
               res_i_g)

    cl.enqueue_copy(queue, res_r, res_r_g)
    cl.enqueue_copy(queue, res_i, res_i_g)
    print(res_r[:8])
    print(res_i[:8])
    print(res_r.sum())
    print(res_i.sum())
    return res_r, res_i


if __name__ == '__main__':
    im = imread('im2.bmp')
    t1 = time.time()
    a = sphericalHarmonics(4)#max(im.shape))
    for i, s in enumerate(a):
        a = (s*im).sum()
        if i < 10:
            print(a)
    print("normal: ", time.time() - t1)
    t2 = time.time()
    res_r, res_i = spharmOCL(im)
    print("OCL: ", time.time() - t2)
    import pdb; pdb.set_trace()
