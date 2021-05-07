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
        Entry(int index, double value): index(index), value(value) { }
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


    // Iterator returns an EntryRef on deference. This is a proxy class for a reference to the
    // underlying value. Value type is Entry.
    class Iterator: public std::iterator<std::random_access_iterator_tag, Entry, int, EntryRef, EntryRef> {
        int *pIndex;
        double *pValue;
    public:
        Iterator(int *indexPtr, double *valuePtr): pIndex(indexPtr), pValue(valuePtr) { }

        Iterator &operator ++() {
            ++pIndex; ++pValue;
            return *this;
        }

        Iterator &operator --() {
            --pIndex; --pValue;
            return *this;
        }

        Iterator operator ++(int) { return Iterator(pIndex++, pValue++); }
        Iterator operator --(int) { return Iterator(pIndex--, pValue--); }

        Iterator &operator +=(int n) {
            pIndex += n; pValue += n;
            return *this;
        }

        Iterator &operator -=(int n) {
            pIndex -= n; pValue -= n;
            return *this;
        }

        bool operator ==(const Iterator &other) const { return pIndex == other.pIndex; }
        bool operator !=(const Iterator &other) const { return pIndex != other.pIndex; }

        Iterator operator +(int n) const { return Iterator(pIndex + n, pValue + n); }
        Iterator operator -(int n) const { return Iterator(pIndex - n, pValue - n); }

        int operator -(const Iterator &other) { return pIndex - other.pIndex; }

        bool operator <(const Iterator &other) { return pIndex < other.pIndex; }

        Iterator operator[](int n) const { return Iterator(pIndex + n, pValue + n); }

        EntryRef operator *() { return EntryRef(*pIndex,*pValue); }
        EntryRef operator ->() { return EntryRef(*pIndex,*pValue); }

    };

    class ConstIterator: public std::iterator<std::random_access_iterator_tag, Entry, int, Entry *, Entry &> {
        int *pIndex;
        double *pValue;
    public:
        ConstIterator(int *indexPtr, double *valuePtr): pIndex(indexPtr), pValue(valuePtr) { }

        ConstIterator &operator ++() {
            ++pIndex; ++pValue;
            return *this;
        }

        ConstIterator &operator --() {
            --pIndex; --pValue;
            return *this;
        }

        ConstIterator operator ++(int) { return ConstIterator(pIndex++, pValue++); }
        ConstIterator operator --(int) { return ConstIterator(pIndex--, pValue--); }

        ConstIterator &operator +=(int n) {
            pIndex += n; pValue += n;
            return *this;
        }

        ConstIterator &operator -=(int n) {
            pIndex -= n; pValue -= n;
            return *this;
        }

        bool operator ==(const ConstIterator &other) const { return pIndex == other.pIndex; }
        bool operator !=(const ConstIterator &other) const { return pIndex != other.pIndex; }

        ConstIterator operator +(int n) const { return ConstIterator(pIndex + n, pValue + n); }
        ConstIterator operator -(int n) const { return ConstIterator(pIndex - n, pValue - n); }

        int operator -(const ConstIterator &other) { return pIndex - other.pIndex; }

        bool operator <(const ConstIterator &other) { return pIndex < other.pIndex; }

        ConstIterator operator[](int n) const { return ConstIterator(pIndex + n, pValue + n); }

        Entry operator *() { return Entry(*pIndex,*pValue); }
        Entry operator ->() { return Entry(*pIndex,*pValue); }
    };



    void test() {
        Iterator myIt = begin();
        Iterator myIt2 = end();
        if(myIt != myIt2) {
            (*myIt).index = 4;
        }

    }

    Iterator begin() { return Iterator(indices.data(), values.data()); }
    Iterator end() { return Iterator(indices.data()+indices.size(), NULL); }

    ConstIterator cbegin() { return ConstIterator(indices.data(), values.data()); }
    ConstIterator cend() { return ConstIterator(indices.data()+indices.size(), NULL); }

protected:
    void swap(SparseVec &rvalue) {
//        std::cout << "Moving " << rvalue.indices << std::endl;
        indices.swap(rvalue.indices);
        values.swap(rvalue.values);
    }

};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);


#endif //GLPKTEST_SPARSEVEC_H
