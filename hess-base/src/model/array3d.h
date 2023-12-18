#pragma once

#include <vector>
using std::vector;

template <typename T>
class Array3D{
public:
    int size_x = 0, size_y = 0, size_z = 0;
    vector<T> data;
    Array3D(int x=0, int y=0, int z=0): size_x(x), size_y(y), size_z(z) {
        data.resize(size_x * size_y * size_z);
    }
    ~Array3D(){
        data.clear();
    }
    void resize(int x, int y, int z){
        size_x = x;
        size_y = y;
        size_z = z;
        data.resize(size_x * size_y * size_z);
    }
    int dim(int i) const{
        if (i == 0)
            return size_x;
        if (i == 1)
            return size_y;
        if (i == 2)
            return size_z;
        return 0;
    }
    T& operator()(int i, int j, int k) { return data[i + size_x * (j + size_y * k)]; }
    const T& operator()(int i, int j, int k) const { return data[i + size_x * (j + size_y * k)]; }
};
