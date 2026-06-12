#include <immintrin.h>

#define __BITMASK5 0x1F //0b00011111

void intrinsic_grad(double *q, int ci, 
                    int nnghbrs, int* ckn, double *cf,  // internal cell values
                    double **qg,                        // Pointer to list of pointer to boundary values
                    int nbf,               double *gcf, // boundary cell values
                    double *gradq){

    // const __mmask8 bitmask5 = 0b00011111; 
    const __m512i vind_x = _mm512_set_epi64 ( 0, 0, 0, 12,  9, 6, 3, 0);
    const __m512i vind_y = _mm512_set_epi64 ( 0, 0, 0, 13, 10, 7, 4, 1);
    const __m512i vind_z = _mm512_set_epi64 ( 0, 0, 0, 14, 11, 8, 5, 2);

    __m512d qi = _mm512_maskz_loadu_pd (__BITMASK5, &q[(ci-1)*5]);
    __m512d grdQx = {};
    __m512d grdQy = {};
    __m512d grdQz = {};
    
    __m512d qk;
    __m512d dq;

    int ck;

    // Internal cell loop
    for (int kcell =0; kcell < nnghbrs; kcell++){
        ck = ckn[kcell];
        qk = _mm512_maskz_loadu_pd (__BITMASK5, &q[(ck-1)*5]);
        dq = _mm512_sub_pd (qk, qi);   // dq = qk - qi

        grdQx = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 0]) , grdQx);
        grdQy = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 1]) , grdQy);
        grdQz = _mm512_fmadd_pd(dq, _mm512_set1_pd (cf[(kcell*3) + 2]) , grdQz);

    }

    for (int kcell =0; kcell < nbf; kcell++){
        qk = _mm512_maskz_loadu_pd (__BITMASK5, qg[kcell]);
        dq = _mm512_sub_pd (qk, qi);   // dq = qk - qi

        grdQx = _mm512_fmadd_pd(dq, _mm512_set1_pd (gcf[(kcell*3) + 0]) , grdQx);
        grdQy = _mm512_fmadd_pd(dq, _mm512_set1_pd (gcf[(kcell*3) + 1]) , grdQy);
        grdQz = _mm512_fmadd_pd(dq, _mm512_set1_pd (gcf[(kcell*3) + 2]) , grdQz);

    }

    // Transpose (using scatter) and write the gradient back to memory
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_x, grdQx, 8);
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_y, grdQy, 8);
    _mm512_mask_i64scatter_pd (&gradq[0], __BITMASK5, vind_z, grdQz, 8);
    
}
