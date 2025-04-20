// waterlib.c
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>

// fast minimum‑image convention
static inline void min_image(double *dx, double *dy, double *dz, const double *box) {
    if (*dx >  0.5*box[0]) *dx -= box[0];
    else if(*dx < -0.5*box[0]) *dx += box[0];
    if (*dy >  0.5*box[1]) *dy -= box[1];
    else if(*dy < -0.5*box[1]) *dy += box[1];
    if (*dz >  0.5*box[2]) *dz -= box[2];
    else if(*dz < -0.5*box[2]) *dz += box[2];
}

static PyObject* py_triplet_angles(PyObject *self, PyObject *args) {
    PyObject *subObj, *posObj, *boxObj;
    double lowCut, highCut;
    if (!PyArg_ParseTuple(args, "OOOdd", &subObj, &posObj, &boxObj, &lowCut, &highCut))
        return NULL;

    // Convert to C-contiguous double arrays
    PyArrayObject *subArr = (PyArrayObject*)PyArray_FROM_OTF(subObj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *posArr = (PyArrayObject*)PyArray_FROM_OTF(posObj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *boxArr = (PyArrayObject*)PyArray_FROM_OTF(boxObj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!subArr || !posArr || !boxArr) {
        Py_XDECREF(subArr); Py_XDECREF(posArr); Py_XDECREF(boxArr);
        PyErr_SetString(PyExc_ValueError,"Could not convert inputs to double arrays");
        return NULL;
    }

    // Check dimensions
    if (PyArray_NDIM(subArr)!=2 || PyArray_DIM(subArr,1)!=3 ||
        PyArray_NDIM(posArr)!=2 || PyArray_DIM(posArr,1)!=3 ||
        PyArray_NDIM(boxArr)!=1 || PyArray_DIM(boxArr,0)!=3)
    {
        Py_DECREF(subArr); Py_DECREF(posArr); Py_DECREF(boxArr);
        PyErr_SetString(PyExc_ValueError,
            "Expected subPos(N,3), Pos(M,3), BoxDims(3,)");
        return NULL;
    }

    npy_intp Ns = PyArray_DIM(subArr,0);
    npy_intp Np = PyArray_DIM(posArr,0);
    double *sub = (double*)PyArray_DATA(subArr);
    double *pos = (double*)PyArray_DATA(posArr);
    double *box = (double*)PyArray_DATA(boxArr);

    double low2  = lowCut*lowCut;
    double high2 = highCut*highCut;

    // 1) First pass: count number of angles
    npy_intp *counts = malloc(sizeof(npy_intp)*Ns);
    if (!counts) { PyErr_NoMemory(); goto fail1; }
    npy_intp total = 0;
    for (npy_intp i=0; i<Ns; ++i) {
        double *pi = sub + 3*i;
        npy_intp neigh=0;
        for (npy_intp j=0; j<Np; ++j) {
            double dx=pos[3*j]-pi[0], dy=pos[3*j+1]-pi[1], dz=pos[3*j+2]-pi[2];
            min_image(&dx,&dy,&dz,box);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > low2 && r2 <= high2) ++neigh;
        }
        counts[i] = neigh*(neigh-1)/2;
        total += counts[i];
    }

    // 2) Allocate outputs exactly
    npy_intp dimsA[1] = { total };
    PyArrayObject *angArr = (PyArrayObject*)PyArray_SimpleNew(1, dimsA, NPY_DOUBLE);
    if (!angArr) { PyErr_NoMemory(); goto fail2; }
    npy_intp dimsC[1] = { Ns };
    PyArrayObject *cntArr = (PyArrayObject*)PyArray_SimpleNew(1, dimsC, NPY_INTP);
    if (!cntArr) { Py_DECREF(angArr); PyErr_NoMemory(); goto fail2; }

    double   *outA = (double*)   PyArray_DATA(angArr);
    npy_intp *outC = (npy_intp*) PyArray_DATA(cntArr);

    // 3) Second pass: fill angles
    npy_intp idx = 0;
    for (npy_intp i=0; i<Ns; ++i) {
        double *pi = sub + 3*i;
        // gather neighbor vectors
        npy_intp neigh=0;
        for (npy_intp j=0; j<Np; ++j) {
            double dx=pos[3*j]-pi[0], dy=pos[3*j+1]-pi[1], dz=pos[3*j+2]-pi[2];
            min_image(&dx,&dy,&dz,box);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > low2 && r2 <= high2) ++neigh;
        }
        double *vec = malloc(sizeof(double)*3*neigh);
        if (!vec) { PyErr_NoMemory(); goto fail3; }

        // fill neighbor displacements
        npy_intp k=0;
        for (npy_intp j=0; j<Np; ++j) {
            double dx=pos[3*j]-pi[0], dy=pos[3*j+1]-pi[1], dz=pos[3*j+2]-pi[2];
            min_image(&dx,&dy,&dz,box);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 <= low2 || r2 > high2) continue;
            vec[3*k]=dx; vec[3*k+1]=dy; vec[3*k+2]=dz; ++k;
        }

        // compute angles
        double inv180pi = 180.0/M_PI;
        for (npy_intp a=0; a<neigh; ++a) {
            double *va = vec+3*a;
            double na = sqrt(va[0]*va[0]+va[1]*va[1]+va[2]*va[2]);
            for (npy_intp b=a+1; b<neigh; ++b) {
                double *vb = vec+3*b;
                double nb = sqrt(vb[0]*vb[0]+vb[1]*vb[1]+vb[2]*vb[2]);
                if (na==0.0 || nb==0.0) {
                    outA[idx++] = 0.0;
                    continue;
                }
                double dot = va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
                double ca  = dot/(na*nb);
                if (ca>1.0) ca=1.0;
                else if (ca<-1.0) ca=-1.0;
                outA[idx++] = acos(ca)*inv180pi;
            }
        }
        free(vec);
        outC[i] = counts[i];
    }

    free(counts);
    Py_DECREF(subArr); Py_DECREF(posArr); Py_DECREF(boxArr);
    return Py_BuildValue("NN", angArr, cntArr);

fail3:
    Py_DECREF(angArr);
    Py_DECREF(cntArr);
fail2:
    free(counts);
fail1:
    Py_DECREF(subArr); Py_DECREF(posArr); Py_DECREF(boxArr);
    return NULL;
}

// module definition
static PyMethodDef WLMeth[] = {
    { "triplet_angles", py_triplet_angles, METH_VARARGS,
      "triplet_angles(subPos(N,3), Pos(M,3), BoxDims(3), lowCut, highCut) -> (angles, counts)" },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef WMod = {
    PyModuleDef_HEAD_INIT, "waterlib", NULL, -1, WLMeth
};

PyMODINIT_FUNC PyInit_waterlib(void) {
    PyObject *m = PyModule_Create(&WMod);
    if (!m) return NULL;
    import_array();
    return m;
}

