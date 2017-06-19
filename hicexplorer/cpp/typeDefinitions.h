#ifndef PYTHON_H
#define PYTHON_H
#include <Python.h>
#endif // PYTHON_H


#ifndef TYPE_DEFINTIONS_H
#define TYPE_DEFINTIONS_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

struct matrixElement {
    uint32_t x;
    uint32_t y;
    int32_t data;
};

#endif // TYPE_DEFINTIONS_H