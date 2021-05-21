//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H

class SparseVec {
public:
    std::vector<int> indices;   // one-based array of indices of non-zero elements (zero'th element stores number of non-zeros)
    std::vector<double> values; // one-based array of values of non-zero elements (zero'th element stores dimension of this vector)

    SparseVec() { }

    SparseVec(int sparseSize): indices(sparseSize), values(sparseSize) { }

    SparseVec(SparseVec &&rvalue) { // move semantics
        swap(rvalue);
    }

    SparseVec(const SparseVec &lvalue): indices(lvalue.indices), values(lvalue.values) { // copy semantics
    }

    SparseVec(const std::vector<double> &dense) {
        for(int i=0; i < dense.size(); ++i) {
            if(double v = dense[i]; v != 0.0) {
                indices.push_back(i);
                values.push_back(v);
            }
        }
    }

//    double operator [](int i) const;
//    double &operator [](int i);

    int sparseSize() const { return indices.size(); }

    void toDense(double *denseVec, int size) const;
    void add(int i, double v);
    void clear();
//    int capacity() const { return std::min(indices.capacity(),values.capacity()); }
//    int dimension() const { return (int)values[0]; }
    void resize(size_t size) { indices.resize(size); values.resize(size); }
    void reserve(size_t n) { indices.reserve(n); values.reserve(n);}

    int *glpkIndexArray() { return indices.data() - 1; }
    double *glpkValueArray() { return values.data() - 1; }

    const int *glpkIndexArray() const { return indices.data() - 1; }
    const double *glpkValueArray() const { return values.data() - 1; }

    //    void setDimension(int dim) { values[0] = (double)dim; }
//    double dotProd(double *dense) const;

    int maxNonZeroIndex() const;

    void sort() {
        std::sort(begin(), end(), [](Entry a, Entry b) {
            return a.index < b.index;
        });
    }


    double operator[](int index) {
        for(int i = 0;i<sparseSize(); ++i) {
            if(indices[i] == index) return values[i];
        }
        return 0.0;
    }

    SparseVec &operator +=(const SparseVec &other) {
        for(int i=0; i < other.sparseSize(); ++i) {
            add(other.indices[i], other.values[i]);
        }
        return *this;
    }

    SparseVec &operator -=(const SparseVec &other) {
        for(int i=0; i < other.sparseSize(); ++i) {
            add(other.indices[i], -other.values[i]);
        }
        return *this;
    }

    double operator *(const std::vector<double> &other) const {
        double dotProd = 0.0;
        for(int i=0; i < sparseSize(); ++i) {
            dotProd += values[i] * other[indices[i]];
        }
        return dotProd;
    }

    friend double operator *(const std::vector<double> &lhs, const SparseVec &rhs) {
        return rhs * lhs;
    }

    friend SparseVec operator +(SparseVec lhs, const SparseVec &rhs) {
        lhs += rhs;
        return lhs; // should use move constructor (or elision?)
    }



    // copy semantics
    SparseVec &operator =(const SparseVec &lvalue) {
        indices = lvalue.indices;
        values = lvalue.values;
        return *this;
    }

    // move semantics
    SparseVec &operator =(SparseVec &&rvalue) {
        swap(rvalue);
        return *this;
    }

    class Entry {
    public:
        const int index;
        const double value;
        Entry(const int index, const double value): index(index), value(value) { }
        Entry *operator ->() { return this; }

        friend std::ostream &operator <<(std::ostream &out, const Entry &entry) {
            out << entry.index << "->" << entry.value;
            return out;
        }

    };

    class EntryRef {
    public:
        int &index;
        double &value;

        EntryRef(int &indexPtr, double &valuePtr): index(indexPtr), value(valuePtr) { }
        EntryRef(const EntryRef &other): EntryRef(other.index, other.value) { };
        EntryRef(EntryRef &&other) noexcept : EntryRef(other.index, other.value) { }

        EntryRef &operator =(const Entry &entry) {
            index = entry.index;
            value = entry.value;
            return *this;
        }

        EntryRef &operator =(const EntryRef &other) {
            index = other.index;
            value = other.value;
            return *this;
        }

        EntryRef &operator =(EntryRef &&other) noexcept {
            index = other.index;
            value = other.value;
            return *this;
        }

        EntryRef *operator->() { return this; }

        operator Entry() const { return Entry(index, value); }

        friend void swap(EntryRef a, EntryRef b) {
            Entry tmp = a;
            a = b;
            b = tmp;
        }

        friend std::ostream &operator <<(std::ostream &out, const EntryRef &entryRef) {
            out << entryRef.index << "->" << entryRef.value;
            return out;
        }
    };


    template<typename DERIVED, typename INTPTR, typename DOUBLEPTR>
    class IteratorBase {
    protected:
        INTPTR      pIndex;
        DOUBLEPTR   pValue;
    public:
        typedef std::random_access_iterator_tag     iterator_category;
        typedef  Entry                              value_type;
        typedef  ptrdiff_t                          difference_type;


        IteratorBase(INTPTR indexPtr, DOUBLEPTR valuePtr): pIndex(indexPtr), pValue(valuePtr) { }

        DERIVED &operator ++() {
            ++pIndex; ++pValue;
            return *(DERIVED *)this;
        }

        DERIVED &operator --() {
            --pIndex; --pValue;
            return *(DERIVED *)this;
        }

        DERIVED operator ++(int) { return DERIVED(pIndex++, pValue++); }
        DERIVED operator --(int) { return DERIVED(pIndex--, pValue--); }

        DERIVED &operator +=(int n) {
            pIndex += n; pValue += n;
            return *(DERIVED *)this;
        }

        DERIVED &operator -=(int n) {
            pIndex -= n; pValue -= n;
            return *(DERIVED *)this;
        }

        bool operator ==(const IteratorBase &other) const { return pIndex == other.pIndex; }
        bool operator !=(const IteratorBase &other) const { return pIndex != other.pIndex; }

        DERIVED operator +(int n) const { return DERIVED(pIndex + n, pValue + n); }
        DERIVED operator -(int n) const { return DERIVED(pIndex - n, pValue - n); }
        friend DERIVED operator +(int n, const DERIVED  &it) { return it + n; }
        difference_type operator -(const IteratorBase &other) { return pIndex - other.pIndex; }

        bool operator <(const IteratorBase &other) { return pIndex < other.pIndex; }
        bool operator >(const IteratorBase &other) { return pIndex > other.pIndex; }
        bool operator <=(const IteratorBase &other) { return pIndex <= other.pIndex; }
        bool operator >=(const IteratorBase &other) { return pIndex >= other.pIndex; }
    };

    // Iterator returns an EntryRef on deference. This is a proxy class for a reference to the
    // underlying value. Value type is Entry.
    class Iterator: public IteratorBase<Iterator, int *, double *> {
    public:
        typedef EntryRef    pointer;
        typedef EntryRef    reference;
        Iterator(int *indexPtr, double *valuePtr): IteratorBase(indexPtr,valuePtr) { }

        EntryRef operator[](int n) const { return EntryRef(*(pIndex + n), *(pValue + n)); }
        EntryRef operator *() const { return EntryRef(*pIndex,*pValue); }
        EntryRef operator ->() const { return EntryRef(*pIndex,*pValue); }

    };

    class ConstIterator: public IteratorBase<ConstIterator, const int *, const double *> {
    public:
        typedef Entry *   pointer;
        typedef Entry &   reference;

        ConstIterator(const int *indexPtr, const double *valuePtr): IteratorBase(indexPtr,valuePtr) { }

        Entry operator[](int n) const { return Entry(*(pIndex + n), *(pValue + n)); }
        Entry operator *() const { return Entry(*pIndex,*pValue); }
        Entry operator ->() const { return Entry(*pIndex,*pValue); }
    };

    Iterator begin() { return Iterator(indices.data(), values.data()); }
    Iterator end() { return Iterator(indices.data()+indices.size(), NULL); }

    ConstIterator begin() const { return cbegin(); }
    ConstIterator end() const { return cend(); }


    ConstIterator cbegin() const { return ConstIterator(indices.data(), values.data()); }
    ConstIterator cend() const { return ConstIterator(indices.data()+indices.size(), NULL); }

protected:
    void swap(SparseVec &rvalue) {
//        std::cout << "Moving " << rvalue.indices << std::endl;
        indices.swap(rvalue.indices);
        values.swap(rvalue.values);
    }

};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);


#endif //GLPKTEST_SPARSEVEC_H
