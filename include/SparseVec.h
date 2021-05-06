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

//    double operator [](int i) const;
//    double &operator [](int i);

    int sparseSize() const { return indices.size(); }

    void toDense(double *denseVec, int size) const;
    void add(int i, double v);
    void clear();
//    int capacity() const { return std::min(indices.capacity(),values.capacity()); }
//    int dimension() const { return (int)values[0]; }
    void resize(int size) { indices.resize(size); values.resize(size); }

    int *glpkIndexArray() { return indices.data() - 1; }
    double *glpkValueArray() { return values.data() - 1; }

    const int *glpkIndexArray() const { return indices.data() - 1; }
    const double *glpkValueArray() const { return values.data() - 1; }

    //    void setDimension(int dim) { values[0] = (double)dim; }
//    double dotProd(double *dense) const;

    int maxNonZeroIndex() const;

    void sort() {
        std::sort(begin(), end(), [](Entry a, Entry b) {
            return a.index() < b.index();
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
    protected:
        int _index;
        double _value;

    public:
        Entry(int index, double value): _index(index), _value(value) { }

        int index() const { return _index; }
        double value() const { return _value; }

    };

    class EntryRef {
    protected:
        int *pIndex;
        double *pValue;

    public:
        EntryRef(int *indexPtr, double *valuePtr): pIndex(indexPtr), pValue(valuePtr) { }

        EntryRef(const EntryRef &other) = default;

        EntryRef(EntryRef &&other) noexcept : pIndex(other.pIndex), pValue(other.pValue) {}

        EntryRef &operator =(const Entry &entry) {
            *pIndex = entry.index();
            *pValue = entry.value();
            return *this;
        }

        EntryRef &operator =(const EntryRef &other) {
            *pIndex = other.index();
            *pValue = other.value();
            return *this;
        }

        EntryRef &operator =(EntryRef &&other) noexcept {
            *pIndex = other.index();
            *pValue = other.value();
            return *this;
        }

        operator Entry() const {
            return Entry(*pIndex, *pValue);
        }

        int index() const { return *pIndex; }
        double value() const { return *pValue; }
        double &value() { return *pValue; }
        int &index() { return *pIndex; }
    };


    // Iterator type. Has the unusual character that it returns itself on deference
    // (the index and value can be accessed/modified via index() and value() members).
    // However, the value_type of the iterator is Entry. STL only specifies that the
    // value returned by deference is convertible to value_type, which is satisfied.
    class Iterator: public std::iterator<std::random_access_iterator_tag, Entry, int, EntryRef *, EntryRef>, EntryRef {
//        int *pIndex;
//        double *pValue;
    public:
        Iterator(int *indexPtr, double *valuePtr): EntryRef(indexPtr, valuePtr) { }
        Iterator(const Iterator &other): EntryRef(other.pIndex, other.pValue) { }

//        int index() const { return *pIndex; }
//        double value() const { return *pValue; }
//        int &index() { return *pIndex; }
//        double &value() { return *pValue; }

        Iterator &operator ++() {
            ++pIndex;
            ++pValue;
            return *this;
        }

        Iterator &operator --() {
            --pIndex;
            --pValue;
            return *this;
        }

        Iterator operator ++(int) { return Iterator(pIndex++, pValue++); }
        Iterator operator --(int) { return Iterator(pIndex--, pValue--); }

        Iterator &operator +=(int n) {
            pIndex += n;
            pValue += n;
            return *this;
        }

        Iterator &operator -=(int n) {
            pIndex -= n;
            pValue -= n;
            return *this;
        }

        bool operator ==(const Iterator &other) const { return pIndex == other.pIndex; }
        bool operator !=(const Iterator &other) const { return pIndex != other.pIndex; }

        Iterator operator +(int n) const { return Iterator(pIndex + n, pValue + n); }
        Iterator operator -(int n) const { return Iterator(pIndex - n, pValue - n); }

        int operator -(const Iterator &other) { return pIndex - other.pIndex; }

        bool operator <(const Iterator &other) { return pIndex < other.pIndex; }

        Iterator operator[](int n) const { return Iterator(pIndex + n, pValue + n); }

//        Iterator &operator =(const Entry &entry) {
//            *pIndex = entry.index();
//            *pValue = entry.value();
//            return *this;
//        }

        Iterator &operator =(const Iterator &other) {
            pIndex = other.pIndex;
            pValue = other.pValue;
            return *this;
        }

//        operator Entry() const { return Entry(*pIndex,*pValue); }

//        const Iterator &operator *() const { return *this; }
        const EntryRef &operator *() const { return *this; }
        const EntryRef *operator ->() const { return this; }
        EntryRef &operator *() { return *this; }
        EntryRef *operator ->() { return this; }

    };

    void test() {
        Iterator myIt = begin();
        Iterator myIt2 = end();
        if(myIt != myIt2) {
            (*myIt).index() = 4;
        }
    }

    Iterator begin() { return Iterator(indices.data(), values.data()); }
    Iterator end() { return Iterator(indices.data()+indices.size(), NULL); }

protected:
    void swap(SparseVec &rvalue) {
//        std::cout << "Moving " << rvalue.indices << std::endl;
        indices.swap(rvalue.indices);
        values.swap(rvalue.values);
    }

};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
